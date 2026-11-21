/* O1a probe shim (own work, 2026-11-20): logs and fails write(2)/writev(2) on the files whose
 * path ends with $SHIM_PATH. Used two ways: LD_PRELOAD for makeblastdb, and linked into ofprobe.
 *   SHIM_PATH   suffix of the tracked file names
 *   SHIM_LIMIT  bytes the tracked fd may take in total (default unlimited)
 *   SHIM_MODE   0: past the limit the write fills up to the limit and returns that short count,
 *                  the next calls fail ENOSPC (a device that fills up)
 *               1: a call that would pass the limit fails ENOSPC and writes nothing
 *   SHIM_LOG    file the event log goes to (appended)
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/syscall.h>
#include <sys/uio.h>
#include <unistd.h>

#define MAXFD 4096
static unsigned char tracked[MAXFD];
static unsigned char injected[MAXFD];   /* the fd the limit applies to */
static char names[MAXFD][64];           /* basename of the tracked path */
static long long total[MAXFD];
static int logfd = -2;
static long long g_limit = -1;
static int g_mode = 0;
static int g_have_cfg = 0;
static int g_close_err = 0;   /* errno close(2) fails with on tracked fds (0: no failure) */

static void cfg(void) {
    if (g_have_cfg) return;
    g_have_cfg = 1;
    const char *l = getenv("SHIM_LIMIT");
    if (l && *l) g_limit = atoll(l);
    const char *ce = getenv("SHIM_CLOSE_ERR");
    if (ce && *ce) g_close_err = atoi(ce);
    const char *m = getenv("SHIM_MODE");
    if (m && *m) g_mode = atoi(m);
}
void shim_close_err(int e) { g_close_err = e; }
void shim_set(long long limit, int mode) { g_have_cfg = 1; g_limit = limit; g_mode = mode; }
static void logline(const char *fmt, ...) {
    if (logfd == -2) {
        const char *p = getenv("SHIM_LOG");
        logfd = p ? (int)syscall(SYS_openat, AT_FDCWD, p, O_WRONLY | O_CREAT | O_APPEND, 0644) : -1;
    }
    if (logfd < 0) return;
    char buf[256];
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    syscall(SYS_write, logfd, buf, n);
}
void shim_mark(const char *s) { logline("%s\n", s); }

/* does `path` end with one of the comma-separated suffixes of the variable `var`? */
static int ends_with_any(const char *var, const char *path) {
    const char *list = getenv(var);
    if (!list || !*list || !path) return 0;
    size_t a = strlen(path);
    while (*list) {
        const char *e = strchr(list, ',');
        size_t b = e ? (size_t)(e - list) : strlen(list);
        if (b && a >= b && memcmp(path + a - b, list, b) == 0) return 1;
        if (!e) break;
        list = e + 1;
    }
    return 0;
}
static int matches(const char *path) { return ends_with_any("SHIM_PATH", path); }
static int track(int fd, const char *path, int flags) {
    if (fd >= 0 && fd < MAXFD && (flags & (O_WRONLY | O_RDWR)) && matches(path)) {
        cfg();
        tracked[fd] = 1;
        total[fd] = 0;
        const char *base = strrchr(path, '/');
        base = base ? base + 1 : path;
        snprintf(names[fd], sizeof names[fd], "%s", base);
        injected[fd] = getenv("SHIM_FAIL") ? ends_with_any("SHIM_FAIL", path) : 1;
        logline("%s open\n", names[fd]);
    }
    return fd;
}
int open(const char *path, int flags, ...) {
    mode_t mode = 0;
    if (flags & O_CREAT) { va_list ap; va_start(ap, flags); mode = va_arg(ap, mode_t); va_end(ap); }
    int fd = (int)syscall(SYS_openat, AT_FDCWD, path, flags, mode);
    if (fd < 0) { errno = -fd; if (matches(path)) logline("open %s FAILED\n", path); return -1; }
    return track(fd, path, flags);
}
int open64(const char *path, int flags, ...) {
    mode_t mode = 0;
    if (flags & O_CREAT) { va_list ap; va_start(ap, flags); mode = va_arg(ap, mode_t); va_end(ap); }
    return open(path, flags, mode);
}
int openat(int dfd, const char *path, int flags, ...) {
    mode_t mode = 0;
    if (flags & O_CREAT) { va_list ap; va_start(ap, flags); mode = va_arg(ap, mode_t); va_end(ap); }
    long r = syscall(SYS_openat, dfd, path, flags, mode);
    if (r < 0) { if (matches(path)) logline("open %s FAILED errno=%d\n", path, errno); return -1; }
    return track((int)r, path, flags);
}
int openat64(int dfd, const char *path, int flags, ...) {
    mode_t mode = 0;
    if (flags & O_CREAT) { va_list ap; va_start(ap, flags); mode = va_arg(ap, mode_t); va_end(ap); }
    return openat(dfd, path, flags, mode);
}
int close(int fd) {
    if (fd >= 0 && fd < MAXFD && tracked[fd]) {
        tracked[fd] = 0;
        int inj = injected[fd] && g_close_err;
        if (inj) logline("%s close -> -1 E%d\n", names[fd], g_close_err); else logline("%s close\n", names[fd]);
        long r = syscall(SYS_close, fd);   /* the descriptor is released either way (Linux) */
        if (inj) { errno = g_close_err; return -1; }
        return (int)r;
    }
    return (int)syscall(SYS_close, fd);
}
/* how many bytes of a call of `count` bytes the shim lets through, or -1 for a refusal */
static long allow(int fd, long long count) {
    if (g_limit < 0 || !injected[fd]) return count;
    long long room = g_limit - total[fd];
    if (room >= count) return count;
    if (g_mode == 0 && room > 0) return room;
    return -1;
}
ssize_t write(int fd, const void *buf, size_t count) {
    if (fd < 0 || fd >= MAXFD || !tracked[fd]) return syscall(SYS_write, fd, buf, count);
    long a = allow(fd, (long long)count);
    if (a < 0) { logline("%s w %zu -> -1 ENOSPC\n", names[fd], count); errno = ENOSPC; return -1; }
    long r = syscall(SYS_write, fd, buf, (size_t)a);
    if (r > 0) total[fd] += r;
    logline("%s w %zu -> %ld\n", names[fd], count, r);
    if (r < 0) { return -1; }
    return r;
}
ssize_t writev(int fd, const struct iovec *iov, int iovcnt) {
    if (fd < 0 || fd >= MAXFD || !tracked[fd]) return syscall(SYS_writev, fd, iov, iovcnt);
    long long count = 0;
    for (int i = 0; i < iovcnt; i++) count += iov[i].iov_len;
    long a = allow(fd, count);
    if (a < 0) { logline("%s v %lld -> -1 ENOSPC\n", names[fd], count); errno = ENOSPC; return -1; }
    if (a < count) {
        /* short: emit the first `a` bytes through plain writes */
        long done = 0;
        for (int i = 0; i < iovcnt && done < a; i++) {
            long take = (long)iov[i].iov_len; if (take > a - done) take = a - done;
            long r = syscall(SYS_write, fd, iov[i].iov_base, (size_t)take);
            if (r < 0) break;
            done += r;
        }
        total[fd] += done;
        logline("%s v %lld -> %ld\n", names[fd], count, done);
        return done;
    }
    long r = syscall(SYS_writev, fd, iov, iovcnt);
    if (r > 0) total[fd] += r;
    logline("%s v %lld -> %ld\n", names[fd], count, r);
    return r;
}
off_t lseek(int fd, off_t off, int whence) {
    long r = syscall(SYS_lseek, fd, off, whence);
    if (fd >= 0 && fd < MAXFD && tracked[fd]) logline("%s l %lld %d -> %ld\n", names[fd], (long long)off, whence, r);
    return r;
}
off_t lseek64(int fd, off_t off, int whence) { return lseek(fd, off, whence); }
ssize_t pwrite(int fd, const void *buf, size_t count, off_t off) {
    if (fd >= 0 && fd < MAXFD && tracked[fd]) logline("%s pwrite %zu\n", names[fd], count);
    return syscall(SYS_pwrite64, fd, buf, count, off);
}

/* the C library's stream open and close (an ofstream opens its file through fopen) */
FILE *fopen(const char *path, const char *mode) {
    FILE *(*real)(const char *, const char *) = (FILE *(*)(const char *, const char *))dlsym(RTLD_NEXT, "fopen");
    FILE *fp = real(path, mode);
    if (!fp) { if (matches(path)) logline("open %s FAILED errno=%d\n", path, errno); return NULL; }
    if (strchr(mode, 'w') || strchr(mode, 'a') || strchr(mode, '+')) track(fileno(fp), path, O_WRONLY);
    return fp;
}
FILE *fopen64(const char *path, const char *mode) { return fopen(path, mode); }
int fclose(FILE *fp) {
    int fd = fileno(fp);
    int inj = 0;
    if (fd >= 0 && fd < MAXFD && tracked[fd]) {
        tracked[fd] = 0;
        inj = injected[fd] && g_close_err;
        if (inj) logline("%s close -> -1 E%d\n", names[fd], g_close_err); else logline("%s close\n", names[fd]);
    }
    int (*real)(FILE *) = (int (*)(FILE *))dlsym(RTLD_NEXT, "fclose");
    int r = real(fp);
    if (inj) { errno = g_close_err; return EOF; }
    return r;
}
