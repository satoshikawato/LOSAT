/* x87 oracle for LOSAT utils/x87.rs: real x87 (gcc long double on x86-64 Linux, control word 0x037F:
 * precision 64 bits, round to nearest even). Reads nothing; writes one line per case:
 *   M e1 s1 e2 s2 e s      fmul        (s = 64-bit significand hex, e = unbiased exponent; zero as "Z")
 *   D e1 s1 e2 s2 e s      fdiv
 *   A hi lo e s            (long double)hi + (long double)lo   (hi, lo = double bits hex)
 *   F e1 s1 bits           fstpl (store to double)
 * usage: x87_oracle <count per op> <seed>                                                           */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
typedef union { long double x; struct { uint64_t m; uint16_t se; } p; } U;
static long double mk(int e, uint64_t sig) { U u; memset(&u, 0, sizeof u); u.p.m = sig; u.p.se = (uint16_t)(e + 16383); return u.x; }
static void put(long double x) { U u; memset(&u, 0, sizeof u); u.x = x;
  if ((u.p.se & 0x7fff) == 0 && u.p.m == 0) { printf(" Z"); return; }
  if ((u.p.se & 0x7fff) == 0) { printf(" DENORMAL"); return; }
  printf(" %d %016llx", (int)(u.p.se & 0x7fff) - 16383, (unsigned long long)u.p.m); }
static uint64_t s = 88172645463325252ull;
static uint64_t rnd(void) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
static uint64_t sig(void) { /* significand with bit 63 set, mixed shapes */
  uint64_t r = rnd() | (1ull << 63); switch (rnd() % 6) {
  case 0: return r;                                   /* uniform */
  case 1: { int tz = rnd() % 64; return (r >> tz << tz) | (1ull << 63); }   /* trailing zeros: exact/tie products */
  case 2: return ~0ull ^ ((rnd() % 4) << (rnd() % 8)); /* near all ones: carries */
  case 3: return (1ull << 63) | (rnd() % 16);          /* near a power of two */
  case 4: { int k = rnd() % 64; return (r | ((1ull << k) - 1)) | (1ull << 63); } /* trailing ones */
  default: return (r & ~0x7ffull) | (rnd() & 1 ? 0x400 : 0x3ff + (rnd() % 3)); /* 53-bit midpoints */ } }
static double dbl_normal(void) { uint64_t b = (rnd() >> 12) | ((uint64_t)(1 + rnd() % 2046) << 52); double d; memcpy(&d, &b, 8); return d; }
int main(int argc, char **argv) {
  long n = atol(argv[1]); s ^= strtoull(argv[2], 0, 10);
  for (long i = 0; i < n; i++) {
    int e1 = (int)(rnd() % 400) - 200, e2 = (int)(rnd() % 400) - 200; uint64_t s1 = sig(), s2 = sig();
    volatile long double a = mk(e1, s1), b = mk(e2, s2), r;
    r = a * b; printf("M %d %016llx %d %016llx", e1, (unsigned long long)s1, e2, (unsigned long long)s2); put(r); printf("\n");
    r = a / b; printf("D %d %016llx %d %016llx", e1, (unsigned long long)s1, e2, (unsigned long long)s2); put(r); printf("\n");
  }
  for (long i = 0; i < n; i++) {
    double hi = (rnd() % 50 == 0) ? ldexp(dbl_normal(), -1060 + (int)(rnd() % 40)) : dbl_normal(); /* some subnormal hi */
    if (hi < 0) hi = -hi;
    int e; frexp(hi, &e); /* hi in [2^(e-1), 2^e): ulp(hi) = 2^(e-53) for normal hi */
    double lo;
    switch (rnd() % 4) {
    case 0: lo = ldexp((double)(int64_t)(rnd() >> 11) / 9007199254740992.0, e - 54); break;        /* |lo| < ulp/2, random */
    case 1: lo = ldexp(1.0, e - 65 - (int)(rnd() % 2)) * ((rnd() & 1) ? 1 : -1); break;           /* x87 ties and quarter ties */
    case 2: lo = ldexp(1.0 + ldexp((double)(rnd() % 1024), -60), e - 65) * ((rnd() & 1) ? 1 : -1); break; /* just past a tie */
    default: lo = ldexp((double)(rnd() % 4096), e - 140 - (int)(rnd() % 80)) * ((rnd() & 1) ? 1 : -1); break; /* far below */
    }
    if (hi < 2.2250738585072014e-308) lo = 0.0; /* a subnormal hi has no lo in a normalized double-double */
    if (fabs(lo) > ldexp(1.0, e - 54)) lo = 0.0;
    volatile long double h = hi, l = lo, r = h + l; uint64_t hb, lb; memcpy(&hb, &hi, 8); memcpy(&lb, &lo, 8);
    printf("A %016llx %016llx", (unsigned long long)hb, (unsigned long long)lb); put(r); printf("\n");
  }
  for (long i = 0; i < n; i++) {
    int e1; uint64_t s1 = sig(); switch (rnd() % 3) { case 0: e1 = (int)(rnd() % 2046) - 1022; break;
      case 1: e1 = -1022 - (int)(rnd() % 60); break; default: e1 = (rnd() & 1) ? 1023 : -1022 - (int)(rnd() % 4); }
    volatile long double a = mk(e1, s1); volatile double d = (double)a; uint64_t db; memcpy(&db, (const void *)&d, 8);
    printf("F %d %016llx %016llx\n", e1, (unsigned long long)s1, (unsigned long long)db);
  }
  return 0;
}
