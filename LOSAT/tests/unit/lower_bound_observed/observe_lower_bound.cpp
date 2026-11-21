// Observation program for line LIB-S item S4 (2026-11-20, ad-LIB-S-1).
//
// It calls std::lower_bound(first, last, key, comp) with a comparator that records every call it
// receives: the index of the element it was given and the answer it returned. Nothing here reads
// or reproduces the library's code; the program only looks at the calls its own comparator gets.
//
// Modes:
//   tree <n> <start> <total>
//       Enumerates every answer sequence the comparator can give for the sub-range
//       [start, start + n) of a vector of `total` elements, and prints the decision tree that the
//       calls form, in preorder: "p<index>" for a call (then the subtree after answer false, then the
//       subtree after answer true), "r<position>" for a returned position (relative to the vector's
//       begin). A last line "stats" gives counts and checks.
//   data <seed> <cases> <maxlen> [direct]
//       Random collections (sorted, sorted with ties, unsorted, and collections built the way
//       NCBI's CRangeCollection<TSeqPos> builds them with x_CombineWith / x_Subtract, including
//       ranges whose to_open wrapped to 0) with the real comparator of NCBI's PRangeLessPos
//       (R.GetToOpen() <= Pos). One line per lower_bound call:
//       "c <case> <shape> <start> <end> <key> <index>:<answer> ... = <returned position>
//        part=<1 if the answers over the window are all true before all false> first_false=<index>".
//       With "direct", only the random sorted / ties / unsorted shapes.
//   one <start> <key> <to_open...>
//       One call on the given collection, printed as a data line (shape "one").
//
// Build: g++-13 -O2 -std=c++17 -o observe_lower_bound observe_lower_bound.cpp
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

struct Range {
    uint32_t from;
    uint32_t to_open;
};

struct Call {
    long index;
    bool answer;
};

struct Context {
    const Range *base = nullptr;
    std::vector<Call> calls;
    const std::vector<char> *script = nullptr;  // forced answers (tree mode)
    bool scripted = false;
};

struct RecordingLess {
    Context *ctx;
    bool operator()(const Range &r, uint32_t key) const {
        long index = &r - ctx->base;
        bool answer;
        if (ctx->scripted) {
            size_t k = ctx->calls.size();
            answer = k < ctx->script->size() ? (*ctx->script)[k] != 0 : false;
        } else {
            answer = r.to_open <= key;  // NCBI PRangeLessPos: R.GetToOpen() <= Pos
        }
        ctx->calls.push_back({index, answer});
        return answer;
    }
};

static long run_scripted(std::vector<Range> &v, long start, long end, const std::vector<char> &script,
                         Context &ctx) {
    ctx.base = v.data();
    ctx.calls.clear();
    ctx.script = &script;
    ctx.scripted = true;
    RecordingLess less{&ctx};
    auto it = std::lower_bound(v.begin() + start, v.begin() + end, uint32_t(0), less);
    return it - v.begin();
}

struct TreeStats {
    long internal = 0;
    long leaves = 0;
    long max_depth = 0;
    long bad_index = 0;      // a call outside [start, end)
    long repeated_index = 0; // a path that gives the comparator the same element twice
    long nondeterministic = 0;
};

static void tree_node(std::vector<Range> &v, long start, long end, std::vector<char> &prefix, Context &ctx,
                      std::string &out, TreeStats &st) {
    long result = run_scripted(v, start, end, prefix, ctx);
    std::vector<Call> calls = ctx.calls;
    // Determinism: the calls up to the prefix must answer as scripted (they do by construction);
    // the path is the same each time this prefix is replayed.
    {
        Context again;
        long result2 = run_scripted(v, start, end, prefix, again);
        if (result2 != result || again.calls.size() != calls.size()) st.nondeterministic++;
    }
    for (size_t i = 0; i < calls.size(); ++i) {
        if (calls[i].index < start || calls[i].index >= end) st.bad_index++;
        for (size_t j = 0; j < i; ++j)
            if (calls[j].index == calls[i].index) st.repeated_index++;
    }
    if (calls.size() == prefix.size()) {
        out += "r" + std::to_string(result) + " ";
        st.leaves++;
        if ((long)calls.size() > st.max_depth) st.max_depth = (long)calls.size();
        return;
    }
    out += "p" + std::to_string(calls[prefix.size()].index) + " ";
    st.internal++;
    prefix.push_back(0);
    tree_node(v, start, end, prefix, ctx, out, st);
    prefix.back() = 1;
    tree_node(v, start, end, prefix, ctx, out, st);
    prefix.pop_back();
}

static int mode_tree(long n, long start, long total) {
    if (start + n > total) {
        fprintf(stderr, "start + n > total\n");
        return 2;
    }
    std::vector<Range> v(total, Range{0, 0});
    std::vector<char> prefix;
    Context ctx;
    std::string out;
    TreeStats st;
    tree_node(v, start, start + n, prefix, ctx, out, st);
    if (!out.empty()) out.pop_back();
    printf("tree %ld %ld %ld %s\n", n, start, total, out.c_str());
    printf("stats %ld %ld %ld internal=%ld leaves=%ld max_depth=%ld bad_index=%ld repeated_index=%ld "
           "nondeterministic=%ld\n",
           n, start, total, st.internal, st.leaves, st.max_depth, st.bad_index, st.repeated_index,
           st.nondeterministic);
    return 0;
}

// ---- data mode -----------------------------------------------------------------------------

static FILE *g_out;
static long g_case = 0;
static const char *g_shape = "";

static long recorded_lower_bound(std::vector<Range> &v, long start, uint32_t key) {
    Context ctx;
    ctx.base = v.data();
    ctx.scripted = false;
    RecordingLess less{&ctx};
    auto it = std::lower_bound(v.begin() + start, v.end(), key, less);
    long result = it - v.begin();
    fprintf(g_out, "c %ld %s %ld %ld %u", g_case, g_shape, start, (long)v.size(), key);
    for (const Call &c : ctx.calls) fprintf(g_out, " %ld:%d", c.index, c.answer ? 1 : 0);
    // Facts about the input (not about the library): whether the comparator's answers over the
    // whole window are all true before all false (partitioned), and the first element of the
    // window for which it answers false (end if none).
    long first_false = (long)v.size();
    bool partitioned = true;
    for (long i = start; i < (long)v.size(); ++i) {
        bool less_i = v[i].to_open <= key;
        if (!less_i && first_false == (long)v.size()) first_false = i;
        if (less_i && first_false != (long)v.size()) partitioned = false;
    }
    fprintf(g_out, " = %ld part=%d first_false=%ld\n", result, partitioned ? 1 : 0, first_false);
    g_case++;
    return result;
}

// NCBI c++/include/util/range_coll.hpp x_CombineWith / x_Subtract (598d8ae6, public domain),
// with the recording comparator in place of PRangeLessPos.
static void combine_with(std::vector<Range> &v, Range r) {
    uint32_t pos_from = r.from;
    uint32_t pos_to_open = r.to_open;
    if (pos_from == 0) pos_from++;
    long b = recorded_lower_bound(v, 0, pos_from - 1);
    if (b != (long)v.size() && v[b].from <= pos_to_open) {
        long e = recorded_lower_bound(v, b, pos_to_open);
        // CombineWith
        bool r_empty = r.to_open <= r.from;
        bool b_empty = v[b].to_open <= v[b].from;
        if (!r_empty) {
            if (!b_empty) {
                v[b].from = std::min(v[b].from, r.from);
                v[b].to_open = std::max(v[b].to_open, r.to_open);
            } else {
                v[b] = r;
            }
        }
        if (e != (long)v.size() && v[e].from <= pos_to_open) {
            v[b].to_open = v[e].to_open;
            ++e;
        }
        v.erase(v.begin() + b + 1, v.begin() + e);
    } else {
        v.insert(v.begin() + b, r);
    }
}

static void subtract(std::vector<Range> &v, Range r) {
    uint32_t pos_from = r.from;
    uint32_t pos_to_open = r.to_open;
    long b = recorded_lower_bound(v, 0, pos_from);
    if (b == (long)v.size()) return;
    if (v[b].from < pos_from && v[b].to_open > pos_to_open) {
        Range copy = v[b];
        v.insert(v.begin() + b, copy);
        v[b].to_open = pos_from;
        v[b + 1].from = pos_to_open;
    } else {
        if (v[b].from < pos_from) {
            v[b].to_open = pos_from;
            ++b;
        }
        long e = recorded_lower_bound(v, b, pos_to_open);
        if (e != (long)v.size() && v[e].from < pos_to_open) v[e].from = pos_to_open;
        v.erase(v.begin() + b, v.begin() + e);
    }
}

static uint32_t pick_key(std::mt19937_64 &rng, const std::vector<Range> &v) {
    uint64_t k = rng() % 8;
    if (v.empty() || k == 0) return (uint32_t)rng();
    const Range &r = v[rng() % v.size()];
    switch (k) {
        case 1: return r.to_open;
        case 2: return r.to_open - 1;
        case 3: return r.to_open + 1;
        case 4: return r.from;
        case 5: return 0;
        case 6: return 0xFFFFFFFFu;
        default: return r.to_open + (uint32_t)(rng() % 64) - 32;
    }
}

static int mode_data(uint64_t seed, long cases, long maxlen, bool direct_only) {
    std::mt19937_64 rng(seed);
    g_out = stdout;
    long built = 0;
    while (g_case < cases) {
        uint64_t shape = direct_only ? rng() % 4 : rng() % 6;
        long len = (long)(rng() % (uint64_t)(maxlen + 1));
        std::vector<Range> v;
        if (shape <= 3) {
            v.resize(len);
            uint32_t span = (rng() % 2) ? 1000u : 0xFFFFFFFFu;
            for (auto &r : v) {
                r.to_open = (uint32_t)(rng() % ((uint64_t)span + 1));
                r.from = r.to_open - (uint32_t)(rng() % 50);
            }
            if (shape == 0) {
                g_shape = "sorted";
                std::sort(v.begin(), v.end(), [](const Range &a, const Range &b) { return a.to_open < b.to_open; });
            } else if (shape == 1) {
                g_shape = "ties";
                for (auto &r : v) r.to_open = (uint32_t)(rng() % 8);
                std::sort(v.begin(), v.end(), [](const Range &a, const Range &b) { return a.to_open < b.to_open; });
            } else if (shape == 2) {
                g_shape = "unsorted";
            } else {
                g_shape = "unsorted_ties";
                for (auto &r : v) r.to_open = (uint32_t)(rng() % 4);
            }
            long queries = 1 + (long)(rng() % 4);
            for (long q = 0; q < queries && g_case < cases; ++q) {
                long start = v.empty() ? 0 : (long)(rng() % (v.size() + 1));
                if (rng() % 3 == 0) start = 0;
                recorded_lower_bound(v, start, pick_key(rng, v));
            }
        } else {
            // Collections the callers build: x_CombineWith / x_Subtract of random ranges; some
            // ranges come from a query gap: from = to = 2^32 - 1 (to_open wraps to 0) or from or to
            // alone at 2^32 - 1, then s_FixMinusStrandRange.
            g_shape = shape == 4 ? "combine" : "subtract";
            long ops = 1 + len;
            uint32_t span = (rng() % 2) ? 2000u : 200000u;
            for (long i = 0; i < ops && g_case < cases; ++i) {
                uint32_t a = (uint32_t)(rng() % span), b = (uint32_t)(rng() % span);
                uint64_t gap = rng() % 10;
                if (gap == 0) a = 0xFFFFFFFFu;
                if (gap == 1) b = 0xFFFFFFFFu;
                if (gap == 2) a = b = 0xFFFFFFFFu;
                Range r{a, b + 1};  // CRange(from, to): to_open = to + 1 (wraps)
                // s_FixMinusStrandRange: if from > to, Set(to, from)
                if (r.from > (uint32_t)(r.to_open - 1)) r = Range{(uint32_t)(r.to_open - 1), r.from + 1};
                if (shape == 4 || v.empty() || rng() % 2 == 0)
                    combine_with(v, r);
                else
                    subtract(v, r);
            }
            built++;
        }
    }
    fprintf(stderr, "data seed=%llu cases=%ld collections_built=%ld\n", (unsigned long long)seed, g_case, built);
    return 0;
}

// one <start> <key> <to_open...>: one call on the given collection (from = to_open), printed as in
// data mode with shape "one".
static int mode_one(long start, uint32_t key, int count, char **values) {
    std::vector<Range> v;
    for (int i = 0; i < count; ++i) {
        uint32_t t = (uint32_t)strtoul(values[i], nullptr, 10);
        v.push_back(Range{t, t});
    }
    g_out = stdout;
    g_shape = "one";
    recorded_lower_bound(v, start, key);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 4 && strcmp(argv[1], "one") == 0)
        return mode_one(atol(argv[2]), (uint32_t)strtoul(argv[3], nullptr, 10), argc - 4, argv + 4);
    if (argc >= 5 && strcmp(argv[1], "tree") == 0)
        return mode_tree(atol(argv[2]), atol(argv[3]), atol(argv[4]));
    if (argc >= 5 && strcmp(argv[1], "data") == 0)
        return mode_data(strtoull(argv[2], nullptr, 10), atol(argv[3]), atol(argv[4]),
                         argc >= 6 && strcmp(argv[5], "direct") == 0);
    fprintf(stderr, "usage: %s tree <n> <start> <total> | data <seed> <cases> <maxlen> [direct] | one <start> <key> <to_open...>\n", argv[0]);
    return 2;
}
