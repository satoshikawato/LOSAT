// Observation program for O1c (string): prints size()/capacity() of a std::string after each
// operation that LOSAT's model must reproduce (reserve(n), append(n,'c'), push_back('c'),
// clear()), the most heap bytes that were live at once during the operation ("peak"), or the
// what() of the exception thrown. Each operation runs on two strings in step: a plain
// std::string and a std::basic_string with a counting allocator (same member functions); the
// program stops with "MISMATCH" if their size/capacity ever differ.
// Input: one operation list per line, e.g. "R30 A5 A100 P C" (R reserve, A append(n,'c'),
// P push_back('c'), C clear).
// Build: x86_64-conda-linux-gnu-g++ 13.4.0 (conda-forge, build 20) -O2 -std=c++20 -static-libstdc++.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <new>
#include <sstream>
#include <string>
typedef unsigned long long u64;

static u64 g_live = 0, g_peak = 0;
template <class T> struct Counting {
    typedef T value_type;
    Counting() = default;
    template <class U> Counting(const Counting<U> &) {}
    T *allocate(std::size_t n) {
        // The same call std::allocator makes: operator new of n * sizeof(T) bytes.
        T *p = std::allocator<T>().allocate(n);
        g_live += (u64)n * sizeof(T);
        if (g_live > g_peak) g_peak = g_live;
        return p;
    }
    void deallocate(T *p, std::size_t n) {
        g_live -= (u64)n * sizeof(T);
        std::allocator<T>().deallocate(p, n);
    }
    template <class U> bool operator==(const Counting<U> &) const { return true; }
};
typedef std::basic_string<char, std::char_traits<char>, Counting<char>> CString;

int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "maxsize")) {
        std::string s;
        printf("max_size %llu\n", (u64)s.max_size());
        return 0;
    }
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream in(line);
        std::string op, out;
        g_live = g_peak = 0;
        std::string s;
        CString c;
        while (in >> op) {
            char kind = op[0];
            u64 n = op.size() > 1 ? strtoull(op.c_str() + 1, 0, 10) : 0;
            g_peak = g_live;
            std::string what, whatc;
            for (int which = 0; which < 2; ++which) {
                try {
                    if (kind == 'R') { if (which) c.reserve((size_t)n); else s.reserve((size_t)n); }
                    else if (kind == 'A') { if (which) c.append((size_t)n, 'c'); else s.append((size_t)n, 'c'); }
                    else if (kind == 'P') { if (which) c.push_back('c'); else s.push_back('c'); }
                    else if (kind == 'C') { if (which) c.clear(); else s.clear(); }
                } catch (std::exception &e) {
                    (which ? whatc : what) = std::string("!") + e.what();
                }
            }
            if (what != whatc || s.size() != c.size() || s.capacity() != c.capacity()) {
                printf("MISMATCH %s\n", line.c_str());
                return 1;
            }
            char buf[200];
            snprintf(buf, sizeof buf, " %s%s%llu/%llu/%llu", what.c_str(), what.empty() ? "" : "|", (u64)s.size(),
                     (u64)s.capacity(), (u64)g_peak);
            out += buf;
        }
        printf("%s\t%s\n", line.c_str(), out.c_str() + (out.empty() ? 0 : 1));
    }
    return 0;
}
