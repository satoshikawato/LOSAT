// Observation program for O1c (vector): the life of rpsbproc's `x_TMatchRec::flex`
// (a std::vector<size_t>) in CProSite::Match, as the match loop drives it:
//   push_back (to a size), clear(), pop_back(), a copy of the vector pushed onto a
//   std::stack<x_TMatchRec> (copy construction; a std::deque of vectors), and
//   `curr = stack.top(); stack.pop()` (copy assignment, then the element is destroyed).
// Operations (one list per input line):
//   G n  clear(), then push_back until size() == n
//   K    clear()
//   Z    pop_back() (the vector must not be empty)
//   S    push a copy of the vector onto the stack
//   N n  the vector = a new vector of n elements (capacity n; the old buffer is released)
//   O    the vector = a copy of the stack's top, then pop the stack (the stack must not be empty)
// After each operation: "capacity/size/stackcount/stack_records/peak": the vector's capacity and
// size, the number of stack entries, the sum of their capacities (records) and the most heap
// records (8 bytes each) that were live at once during the operation, all buffers of the vector
// and of the stack's copies counted by an allocator that calls the same operator new.
// Build: x86_64-conda-linux-gnu-g++ 13.4.0 -O2 -std=c++20 -static-libstdc++.
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <exception>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
typedef unsigned long long u64;
static u64 g_live = 0, g_peak = 0;
template <class T> struct Counting {
    typedef T value_type;
    Counting() = default;
    template <class U> Counting(const Counting<U> &) {}
    T *allocate(std::size_t n) {
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
typedef std::vector<size_t, Counting<size_t>> Vec;
int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream in(line);
        std::string op, out;
        g_live = g_peak = 0;
        Vec v;
        std::deque<Vec> stack;
        while (in >> op) {
            char kind = op[0];
            u64 n = op.size() > 1 ? strtoull(op.c_str() + 1, 0, 10) : 0;
            g_peak = g_live;
            std::string what;
            try {
                if (kind == 'G') { v.clear(); for (u64 i = 0; i < n; ++i) v.push_back(i); }
                else if (kind == 'K') v.clear();
                else if (kind == 'N') { Vec w(n); v = std::move(w); }
                else if (kind == 'Z') v.pop_back();
                else if (kind == 'S') stack.push_back(v);
                else if (kind == 'O') { v = stack.back(); stack.pop_back(); }
            } catch (std::exception &e) { what = std::string("!") + e.what() + "|"; }
            u64 sc = 0;
            for (auto &e : stack) sc += e.capacity();
            char buf[200];
            snprintf(buf, sizeof buf, " %s%llu/%llu/%llu/%llu/%llu", what.c_str(), (u64)v.capacity(), (u64)v.size(),
                     (u64)stack.size(), sc, g_peak / 8);
            out += buf;
        }
        printf("%s\t%s\n", line.c_str(), out.c_str() + (out.empty() ? 0 : 1));
    }
    return 0;
}
