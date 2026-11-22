// LD_PRELOAD shim for O1c: logs every malloc/calloc/realloc/free of size >= MIN_LOG (env
// MALLOG_MIN, default 16) to the file MALLOG_FILE: "m <size> <ptr>", "r <size> <oldptr> <ptr>",
// "f <size> <ptr>" (free of a pointer that was logged, with its size). Observation only; calls the next allocator.
#define _GNU_SOURCE
#include <dlfcn.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
static void *(*real_malloc)(size_t);
static void *(*real_calloc)(size_t, size_t);
static void *(*real_realloc)(void *, size_t);
static void (*real_free)(void *);
static int fd = -1;
static size_t min_log = 16;
static __thread int busy;
#define N 1048576
static void *seen[N]; // open-addressing set of logged pointers
static size_t sizes[N];
static unsigned long hsh(void *p) { return ((unsigned long)p >> 4) * 11400714819323198485ul >> 44; }
static void put(void *p, size_t n) { unsigned long i = hsh(p) & (N - 1); while (seen[i] && seen[i] != (void *)1) i = (i + 1) & (N - 1); seen[i] = p; sizes[i] = n; }
static long del(void *p) { unsigned long i = hsh(p) & (N - 1); while (seen[i]) { if (seen[i] == p) { seen[i] = (void *)1; return (long)sizes[i]; } i = (i + 1) & (N - 1); } return -1; }
static int initing;
static char arena[1 << 16]; static size_t used;
static void *early(size_t n) { n = (n + 15) & ~15ul; void *p = arena + used; used += n; return p; }
static void init(void) {
    if (real_malloc || initing) return;
    initing = 1;
    real_malloc = dlsym(RTLD_NEXT, "malloc"); real_calloc = dlsym(RTLD_NEXT, "calloc");
    real_realloc = dlsym(RTLD_NEXT, "realloc"); real_free = dlsym(RTLD_NEXT, "free");
    const char *f = getenv("MALLOG_FILE"), *m = getenv("MALLOG_MIN");
    if (m) min_log = strtoull(m, 0, 10);
    if (f) fd = open(f, O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
    initing = 0;
}
static void logline(const char *fmt, size_t a, void *b, void *c) {
    char buf[128]; int n = snprintf(buf, sizeof buf, fmt, a, b, c);
    if (fd >= 0) { ssize_t r = write(fd, buf, n); (void)r; }
}
void *malloc(size_t n) {
    if (!real_malloc) { init(); if (!real_malloc) return early(n); }
    void *p = real_malloc(n);
    if (!busy && n >= min_log) { busy = 1; put(p, n); logline("m %zu %p%.0s\n", n, p, 0); busy = 0; }
    return p;
}
void *calloc(size_t a, size_t b) {
    if (!real_malloc) init();
    if (!real_calloc) return early(a * b);
    void *p = real_calloc(a, b);
    if (!busy && a * b >= min_log) { busy = 1; put(p, a * b); logline("m %zu %p%.0s\n", a * b, p, 0); busy = 0; }
    return p;
}
void *realloc(void *o, size_t n) {
    if (!real_malloc) init();
    if (!real_realloc) return early(n);
    void *p = real_realloc(o, n);
    if (!busy && n >= min_log) { busy = 1; long was = o ? del(o) : -1; put(p, n); logline("r %zu %p %p\n", n, o, p); (void)was; busy = 0; }
    return p;
}
void free(void *p) {
    if ((char *)p >= arena && (char *)p < arena + sizeof arena) return;
    if (!real_malloc) init();
    long sz; if (!busy && p && (sz = del(p)) >= 0) { busy = 1; logline("f %zu %p%.0s\n", (size_t)sz, p, 0); busy = 0; }
    if (real_free) real_free(p);
}
