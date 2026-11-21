// Records every call that std::list<int>::sort makes to a comparator we pass: which two elements
// (by their input positions) it hands over, in which argument order, and the answer we gave, plus
// the order of the list after the sort. Only the program's inputs and outputs are recorded; no
// library source, header text, preprocessed output, disassembly or debugger trace is read.
//
//   tree <n>                     every answer sequence for a list of n elements (a complete
//                                decision tree: replay each prefix, branch each new call on
//                                false and true); one line per leaf
//   rand <seed> <cases> <maxn>   answers drawn at random per call (p(true) varies by case)
//   list <seed> <permille> <n>... one list of each length n, answers true with p = permille/1000
//
// Line formats (fields separated by one space):
//   T <n> <answers as 0/1, '-' if none> <calls a:b,... or '-'> <order a,b,... or '-'>
//   R <case> <n> <p_true_permille> <answers> <calls> <order>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <random>
#include <string>
#include <vector>

struct Run {
    std::vector<int> a, b;   // arguments of each call (input positions)
    std::vector<char> ans;   // answer of each call
    std::vector<int> order;  // list after the sort
};

template <class Answer>
static Run run(int n, Answer answer) {
    Run r;
    std::list<int> l;
    for (int i = 0; i < n; ++i) l.push_back(i);
    auto comp = [&](const int& x, const int& y) {
        bool v = answer(r.ans.size());
        r.a.push_back(x);
        r.b.push_back(y);
        r.ans.push_back(v ? 1 : 0);
        return v;
    };
    l.sort(comp);
    r.order.assign(l.begin(), l.end());
    return r;
}

static void put(const Run& r) {
    std::string s;
    if (r.ans.empty()) s += '-';
    for (char c : r.ans) s += c ? '1' : '0';
    s += ' ';
    if (r.a.empty()) s += '-';
    for (size_t i = 0; i < r.a.size(); ++i) {
        if (i) s += ',';
        s += std::to_string(r.a[i]) + ':' + std::to_string(r.b[i]);
    }
    s += ' ';
    if (r.order.empty()) s += '-';
    for (size_t i = 0; i < r.order.size(); ++i) {
        if (i) s += ',';
        s += std::to_string(r.order[i]);
    }
    std::fputs(s.c_str(), stdout);
    std::fputc('\n', stdout);
}

static long leaves = 0;

// Every path whose answers start with `script`: the run answers false past the script; each call
// past the script is also branched on true.
static void tree(int n, std::vector<char>& script) {
    size_t k = script.size();
    Run r = run(n, [&](size_t i) { return i < script.size() ? script[i] != 0 : false; });
    std::printf("T %d ", n);
    put(r);
    ++leaves;
    for (size_t j = k; j < r.ans.size(); ++j) {
        std::vector<char> next(script);
        next.resize(j, 0);
        next.push_back(1);
        tree(n, next);
    }
}

int main(int argc, char** argv) {
    if (argc >= 3 && std::string(argv[1]) == "tree") {
        int n = std::atoi(argv[2]);
        std::vector<char> script;
        tree(n, script);
        std::fprintf(stderr, "tree n=%d leaves=%ld\n", n, leaves);
        return 0;
    }
    if (argc >= 5 && std::string(argv[1]) == "rand") {
        unsigned long long seed = std::strtoull(argv[2], nullptr, 10);
        long cases = std::atol(argv[3]);
        int maxn = std::atoi(argv[4]);
        std::mt19937_64 rng(seed);
        const int permille[] = {500, 100, 900, 20, 980, 300, 700};
        for (long c = 0; c < cases; ++c) {
            int n = (int)(rng() % (unsigned long long)(maxn + 1));
            if (c % 4 == 0) n = (int)(rng() % 17ULL);  // many short lists
            int p = permille[rng() % 7];
            Run r = run(n, [&](size_t) { return (int)(rng() % 1000ULL) < p; });
            std::printf("R %ld %d %d ", c, n, p);
            put(r);
        }
        return 0;
    }
    if (argc >= 5 && std::string(argv[1]) == "list") {
        std::mt19937_64 rng(std::strtoull(argv[2], nullptr, 10));
        int p = std::atoi(argv[3]);
        for (int i = 4; i < argc; ++i) {
            int n = std::atoi(argv[i]);
            Run r = run(n, [&](size_t) { return (int)(rng() % 1000ULL) < p; });
            std::printf("R %d %d %d ", i - 4, n, p);
            put(r);
        }
        return 0;
    }
    std::fprintf(stderr, "usage: tree <n> | rand <seed> <cases> <maxn> | list <seed> <permille> <n>...\n");
    return 2;
}
