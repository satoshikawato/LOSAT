// O1a probe (own work, 2026-11-20): runs scripted sequences of operations on a std::ofstream and
// prints the stream state after each operation. The syscalls come from shim.c (linked in).
// Script lines (one sequence = R ... Z):
//   R <limit> <mode>   new sequence; the shim's byte limit (-1 none) and mode
//   O                  open <out>             OB  open a path that cannot be opened (a directory)
//   OF open /dev/full  OM  open nothing (stream without a file)
//   W n   write n bytes     I n   operator<< of a std::string of n bytes     P  put one byte
//   F flush   E endl   S pos seekp   C close   X clear(goodbit)   T bits  setstate(bits)
//   Z  end: destroy the stream, print the size and checksum of <out>
// Bytes written come from a running counter: byte k = (k*37 + (k>>8)*11 + 1) & 255.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <dlfcn.h>
// shim_set / shim_mark come from shim.so (LD_PRELOAD); without it they do nothing
static void shim_set(long long l, int m) { auto f = (void (*)(long long, int))dlsym(RTLD_DEFAULT, "shim_set"); if (f) f(l, m); }
static void shim_close_err(int e) { auto f = (void (*)(int))dlsym(RTLD_DEFAULT, "shim_close_err"); if (f) f(e); }
static void shim_mark(const char *s) { auto f = (void (*)(const char *))dlsym(RTLD_DEFAULT, "shim_mark"); if (f) f(s); }
static unsigned long long ctr = 0;
static std::string chunk(size_t n) {
    std::string s(n, '\0');
    for (size_t i = 0; i < n; i++, ctr++) s[i] = (char)((ctr * 37 + (ctr >> 8) * 11 + 1) & 255);
    return s;
}
int main(int argc, char **argv) {
    if (argc < 3) return 2;
    std::ifstream script(argv[1]);
    std::string out = argv[2];
    std::ofstream *f = nullptr;
    std::string line;
    long seq = 0, op = 0;
    while (std::getline(script, line)) {
        if (line.empty()) continue;
        std::istringstream is(line);
        std::string c;
        is >> c;
        char mark[128];
        if (c == "R") {
            long long limit; int mode;
            is >> limit >> mode;
            shim_set(limit, mode);
            ctr = 0; op = 0; seq++; shim_close_err(0);
            snprintf(mark, sizeof mark, "SEQ %ld limit=%lld mode=%d", seq, limit, mode);
            shim_mark(mark);
            f = new std::ofstream();
            continue;
        }
        if (c == "Z") {
            shim_mark("DESTROY");
            delete f; f = nullptr;
            shim_mark("END");
            struct stat st;
            unsigned long long h = 1469598103934665603ULL; long long size = -1;
            if (stat(out.c_str(), &st) == 0) {
                size = st.st_size;
                FILE *g = fopen(out.c_str(), "rb");
                int ch;
                while (g && (ch = fgetc(g)) != EOF) { h ^= (unsigned)ch; h *= 1099511628211ULL; }
                if (g) fclose(g);
            }
            printf("FILE %lld %016llx\n", size, h);
            continue;
        }
        op++;
        snprintf(mark, sizeof mark, "OP %ld %s", op, line.c_str());
        shim_mark(mark);
        if (c == "O") f->open(out.c_str(), std::ios::out | std::ios::binary);
        else if (c == "OB") f->open("/tmp", std::ios::out | std::ios::binary);
        else if (c == "OF") f->open("/dev/full", std::ios::out | std::ios::binary);
        else if (c == "OM") {}
        else if (c == "W") { long n; is >> n; std::string s = chunk(n); f->write(s.data(), (std::streamsize)s.size()); }
        else if (c == "I") { long n; is >> n; std::string s = chunk(n); *f << s; }
        else if (c == "P") { std::string s = chunk(1); f->put(s[0]); }
        else if (c == "F") f->flush();
        else if (c == "E") *f << std::endl;
        else if (c == "S") { long long p; is >> p; f->seekp((std::streamoff)p); }
        else if (c == "C") f->close();
        else if (c == "K") { int e; is >> e; shim_close_err(e); }
        else if (c == "X") f->clear(std::ios::goodbit);
        else if (c == "T") { int b; is >> b; f->setstate((std::ios::iostate)b); }
        printf("STATE %d\n", (int)f->rdstate());
    }
    return 0;
}
