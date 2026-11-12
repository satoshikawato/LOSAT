// DIAGNOSTIC ONLY (branch diag/ph-webkit-hang, line WFL item PH; not for merge). Copied from
// losat-baselines/web-audit-20261104/w1/crash/bt/crashbt.c (line WA W1).
// LD_PRELOAD into WebKit's processes (WPE*, MiniBrowser): each start appends a line to
// $CRASHBT_DIR/starts.log (to tell this run's processes from other sessions' in dmesg), and on a fatal signal, write the signal, the thread name,
// the backtrace addresses and /proc/self/maps to $CRASHBT_DIR/<pid>.<tid>.<sig>.txt, then die
// the same way (fault: default action and re-execute; SIGTRAP/SIGABRT: default and re-raise).
// Other programs (node, the test runner) load it and do nothing.
// With $FAKE_MEMINFO set, MiniBrowser (the UI process, where WebKit's MemoryPressureMonitor
// polls /proc/meminfo) reads that file instead of /proc/meminfo; the web processes read the real one.
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <execinfo.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/syscall.h>
#include <time.h>
#include <ucontext.h>
#include <unistd.h>

static char out_dir[512];
static const int sigs[] = {SIGSEGV, SIGBUS, SIGILL, SIGFPE, SIGTRAP, SIGABRT};

static const char *fake_meminfo;

static FILE *fake_or_real(const char *sym, const char *path, const char *mode) {
  FILE *(*real)(const char *, const char *) = (FILE * (*)(const char *, const char *)) dlsym(RTLD_NEXT, sym);
  if (fake_meminfo && strcmp(path, "/proc/meminfo") == 0) return real(fake_meminfo, mode);
  return real(path, mode);
}
FILE *fopen(const char *path, const char *mode) { return fake_or_real("fopen", path, mode); }
FILE *fopen64(const char *path, const char *mode) { return fake_or_real("fopen64", path, mode); }

static void put(int fd, const char *s) { (void)!write(fd, s, strlen(s)); }

static void put_hex(int fd, const char *label, unsigned long v) {
  char b[64];
  int n = snprintf(b, sizeof b, "%s0x%lx\n", label, v);
  (void)!write(fd, b, n);
}

static void copy_file(int fd, const char *path) {
  int in = open(path, O_RDONLY);
  if (in < 0) return;
  char b[8192];
  ssize_t n;
  while ((n = read(in, b, sizeof b)) > 0) (void)!write(fd, b, n);
  close(in);
}

static void append_line(const char *name, const char *fmt, int a, int b, int c) {
  char path[640], line[256];
  snprintf(path, sizeof path, "%s/%s", out_dir, name);
  int fd = open(path, O_WRONLY | O_CREAT | O_APPEND, 0644);
  if (fd < 0) return;
  int n = snprintf(line, sizeof line, fmt, a, b, c);
  (void)!write(fd, line, n);
  close(fd);
}

static void handler(int sig, siginfo_t *si, void *uc_) {
  ucontext_t *uc = uc_;
  pid_t pid = getpid(), tid = syscall(SYS_gettid);
  char path[640];
  append_line("entered.log", "pid %d tid %d sig %d\n", pid, tid, sig);
  snprintf(path, sizeof path, "%s/%d.%d.%d.txt", out_dir, pid, tid, sig);
  int fd = open(path, O_WRONLY | O_CREAT | O_APPEND, 0644);
  if (fd >= 0) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    char b[256];
    int n = snprintf(b, sizeof b, "time %ld.%09ld\nsignal %d code %d addr %p\npid %d tid %d\n",
                     (long)ts.tv_sec, ts.tv_nsec, sig, si->si_code, si->si_addr, pid, tid);
    (void)!write(fd, b, n);
    snprintf(b, sizeof b, "/proc/self/task/%d/comm", tid);
    put(fd, "comm ");
    copy_file(fd, b);
    put_hex(fd, "rip ", uc->uc_mcontext.gregs[REG_RIP]);
    put_hex(fd, "rsp ", uc->uc_mcontext.gregs[REG_RSP]);
    put_hex(fd, "rdi ", uc->uc_mcontext.gregs[REG_RDI]);
    put_hex(fd, "rsi ", uc->uc_mcontext.gregs[REG_RSI]);
    void *frames[128];
    int depth = backtrace(frames, 128);
    put(fd, "backtrace\n");
    for (int i = 0; i < depth; i++) put_hex(fd, "  ", (unsigned long)frames[i]);
    put(fd, "symbols\n");
    backtrace_symbols_fd(frames, depth, fd);
    put(fd, "status\n");
    copy_file(fd, "/proc/self/status");
    put(fd, "meminfo\n");
    copy_file(fd, "/proc/meminfo");
    put(fd, "maps\n");
    copy_file(fd, "/proc/self/maps");
    close(fd);
  }
  signal(sig, SIG_DFL);
  if (sig == SIGTRAP || sig == SIGABRT || si->si_code <= 0) raise(sig);
}

__attribute__((constructor)) static void init(void) {
  const char *me = program_invocation_short_name;
  if (strncmp(me, "WPE", 3) != 0 && strcmp(me, "MiniBrowser") != 0) return;
  const char *d = getenv("CRASHBT_DIR");
  if (!d) return;
  snprintf(out_dir, sizeof out_dir, "%s", d);
  if (strcmp(me, "MiniBrowser") == 0) fake_meminfo = getenv("FAKE_MEMINFO");
  const char *cpu = getenv("WEBKIT_SKIA_ENABLE_CPU_RENDERING");
  char fmt[160];
  snprintf(fmt, sizeof fmt, "start pid %%d ppid %%d cpu %%d %s\n", me);
  append_line("starts.log", fmt, getpid(), getppid(), cpu ? atoi(cpu) : -1);
  void *warm[4];
  backtrace(warm, 4);  // loads libgcc_s now, not inside the handler
  static char alt[65536];
  stack_t ss = {.ss_sp = alt, .ss_size = sizeof alt};
  sigaltstack(&ss, NULL);
  for (size_t i = 0; i < sizeof sigs / sizeof sigs[0]; i++) {
    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_sigaction = handler;
    sa.sa_flags = SA_SIGINFO | SA_ONSTACK;
    sigemptyset(&sa.sa_mask);
    sigaction(sigs[i], &sa, NULL);
  }
}
