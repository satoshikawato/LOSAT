/* Expected values of the utils/x87.rs unit tests, from the real x87 (gcc long double, x86-64 Linux). */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef union { long double x; struct { uint64_t m; uint16_t se; } p; } U;
static long double mk(int e, uint64_t s) { U u; memset(&u, 0, sizeof u); u.p.m = s; u.p.se = (uint16_t)(e + 16383); return u.x; }
static void put(const char *tag, long double x) { U u; memset(&u, 0, sizeof u); u.x = x;
  printf("%s -> %s%d 0x%016llx\n", tag, (u.p.se & 0x8000) ? "-" : "", (int)(u.p.se & 0x7fff) - 16383, (unsigned long long)u.p.m); }
static void putd(const char *tag, double d) { uint64_t b; memcpy(&b, &d, 8); printf("%s -> 0x%016llx\n", tag, (unsigned long long)b); }
int main(void) {
  volatile long double a, b;
  a = mk(0, 0xC000000000000000ull); b = mk(0, 0x8000000000000001ull); put("mul tie odd", a * b);
  a = mk(0, 0xC000000000000000ull); b = mk(0, 0x8000000000000003ull); put("mul tie even", a * b);
  a = mk(0, 0xFFFFFFFFFFFFFFFEull); b = mk(0, 0x8000000000000001ull); put("mul carry", a * b);
  a = mk(0, 0x8000000000000000ull); b = mk(0, 0xC000000000000000ull); put("div 1/1.5", a / b);
  a = mk(0, 0x8000000000000000ull); b = mk(2, 0xE000000000000000ull); put("div 1/7", a / b);
  a = mk(1, 0xC000000000000000ull); b = mk(0, 0xE000000000000000ull); put("div 3/1.75", a / b);
  a = mk(0, 0xFFFFFFFFFFFFFFFFull); b = mk(0, 0xAAAAAAAAAAAAAAABull); put("div q>=2^64 up?", a / b);
  a = mk(0, 0xFFFFFFFFFFFFFFFFull); b = mk(0, 0xC000000000000001ull); put("div q>=2^64 b", a / b);
  a = -mk(5, 0x9000000000000000ull); b = mk(1, 0xC000000000000000ull); put("div neg 36/3", a / b);
  { volatile long double h, l;
    h = 1.0; l = ldexp(1, -64); put("sum tie even", h + l);
    h = 1.0; l = ldexp(1.5, -63); put("sum tie odd (1+1.5*2^-63+... no)", h + l);
    h = 1.0 + ldexp(1, -52); l = ldexp(1.5, -63); put("sum tie odd", h + l);
    h = 1.0; l = ldexp(1, -64) + ldexp(1, -100); put("sum sticky above tie", h + l);
    h = 1.0; l = -ldexp(1, -65); put("sum borrow tie", h + l);
    h = 1.0; l = -ldexp(1, -64); put("sum borrow exact", h + l);
    h = 1.0; l = -(ldexp(1, -65) + ldexp(1, -117)); put("sum borrow below tie", h + l);
    h = 1.0 + ldexp(1, -52); l = -(ldexp(1, -64) + ldexp(1, -110)); put("sum neg sticky", h + l);
    h = ldexp(1, -1070); l = 0; put("sum subnormal hi", h + l);
  }
  putd("store tie even", (double)mk(0, 0x8000000000000400ull));
  putd("store tie odd", (double)mk(0, 0x8000000000000C00ull));
  putd("store above tie", (double)mk(0, 0x8000000000000401ull));
  putd("store carry", (double)mk(0, 0xFFFFFFFFFFFFFC00ull));
  putd("store overflow", (double)mk(1023, 0xFFFFFFFFFFFFFC00ull));
  putd("store min subnormal", (double)mk(-1074, 0x8000000000000000ull));
  putd("store subnormal tie even", (double)mk(-1075, 0x8000000000000000ull));
  putd("store subnormal tie odd", (double)mk(-1074, 0xC000000000000000ull));
  putd("store subnormal above half", (double)mk(-1075, 0x8000000000000001ull));
  putd("store largest subnormal up", (double)mk(-1023, 0xFFFFFFFFFFFFFFFFull));
  putd("store subnormal 0x...7FF.8", (double)mk(-1023, 0xFFFFFFFFFFFFE800ull));
  putd("store underflow", (double)mk(-1100, 0x8000000000000000ull));
  putd("store negative", (double)(-mk(0, 0xC000000000000400ull)));
  return 0;
}
