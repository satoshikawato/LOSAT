// observe_sort.cpp: the calls that C++ `std::sort` (and `std::partial_sort` over a whole range)
// makes to a comparator, and the order it leaves.
//
// The program sorts the input positions 0..n-1 (ints) with a comparator that writes down every
// call it receives: the position handed as its first argument, the one handed as its second, and
// the answer it gave. After the sort it writes the order of the positions. Nothing else is looked
// at: no library header or source, no preprocessed output, no disassembly or debugger trace.
//
// The n positions sit in the middle of a larger buffer whose other cells hold numbers outside
// 0..n-1. A read outside the sorted range would reach the comparator as such a number; it is
// counted (field `oob`, and a total on standard error) and answered false.
//
// Output, one line per sort:
//   <tag> <alg> <n> <answers> <calls> <order> <oob> [<extra>]
//   tag     T (tree), K (random keys), G (NCBI-shaped ranges), A (lazy adversary), F (keys read
//           from standard input)
//   alg     s = std::sort, p = std::partial_sort(first, last, last); S, P: the same on a heavier
//           element type (Heavy), to show that the calls do not depend on the type
//   answers 0/1 per call, in call order ("-" when there is no call)
//   calls   a:b per call, the positions handed as first and second argument ("-" when none)
//   order   the positions after the sort ("-" when n = 0)
//   oob     number of calls that were handed a cell outside the range
//   extra   K, F: <shape>/<keys>; G: <reverse>/<from:to_open,...>; A: <tie call or -1>
//
// Modes:
//   tree <alg> <n>                     every answer sequence that some strict weak order gives
//   keys <alg> <seed> <cases> <maxn>   random strict weak orders with many ties, several shapes
//   ranges <alg> <seed> <cases> <maxn> NCBI's CRange_Less / CRange_ReverseLess on (from, to_open)
//   adv <alg> <seed> <cases> <maxn>    a comparator that decides lazily (see adversary_answer)
//   file <alg>                         one sort per input line "<n> <key>..." (keys compared with <)
//
// Build: g++-13 -O2 -std=c++17 -o observe_sort observe_sort.cpp
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

const int GUARD = 64;  // cells on each side of the sorted range

enum class Mode { Keys, Tree, Ranges, Adversary };

struct State {
  Mode mode = Mode::Keys;
  int n = 0;
  long oob = 0;
  std::vector<int> calls_a, calls_b;
  std::vector<char> answers;
  // Keys
  std::vector<uint64_t> key;
  // Ranges
  std::vector<uint32_t> from, to_open;
  bool reverse = false;
  // Tree: the answers given so far over all runs (a prefix is replayed), whether the other
  // answer (true) is also consistent at each call, and the closure of what the answers imply.
  std::vector<char> tree_answers, tree_alt;
  size_t tree_prefix = 0;
  std::vector<uint8_t> le, lt;  // n x n: "i <= j", "i < j" implied
  // Adversary
  std::vector<int64_t> value;  // -1 = not decided yet
  int64_t decided = 0;
  int candidate = -1;
  long tie_at = -1;  // call number at which the undecided elements become ties (-1: never)
  std::mt19937_64 tie_rng;
};

State st;

bool range_whole(int i) { return st.from[i] == 0 && st.to_open[i] == 0xFFFFFFFFu; }
bool range_empty(int i) { return st.to_open[i] <= st.from[i]; }

// NCBI's CRange_Less / CRange_ReverseLess (Seq_loc.cpp:4283-4322), on the stored ranges.
bool range_less(int a, int b) {
  if (range_whole(a)) return !range_whole(b);
  if (range_empty(a)) return !range_empty(b) && !range_whole(b);
  if (range_whole(b) || range_empty(b)) return false;
  if (st.reverse) {
    uint32_t ta = st.to_open[a] - 1, tb = st.to_open[b] - 1;
    if (ta != tb) return ta > tb;
    return st.from[a] < st.from[b];
  }
  return st.from[a] < st.from[b] || (st.from[a] == st.from[b] && st.to_open[a] < st.to_open[b]);
}

// Adds "u <= v" (with `strict`: "u < v") to the closure.
void tree_add(int u, int v, bool strict) {
  int n = st.n;
  for (int i = 0; i < n; i++) {
    if (!st.le[i * n + u]) continue;
    bool before = st.lt[i * n + u];
    for (int j = 0; j < n; j++) {
      if (!st.le[v * n + j]) continue;
      st.le[i * n + j] = 1;
      if (strict || before || st.lt[v * n + j]) st.lt[i * n + j] = 1;
    }
  }
}

bool tree_answer(int x, int y) {
  int n = st.n;
  size_t k = st.answers.size();
  bool can_true = !st.le[y * n + x];   // x < y, unless y <= x is implied
  bool can_false = !st.lt[x * n + y];  // y <= x, unless x < y is implied
  bool answer;
  if (k < st.tree_prefix) {
    answer = st.tree_answers[k];
    if ((answer && !can_true) || (!answer && !can_false)) {
      fprintf(stderr, "inconsistent script at call %zu\n", k);
      exit(2);
    }
  } else {
    answer = !can_false;  // false first, when it is consistent
    st.tree_answers.push_back(answer);
    st.tree_alt.push_back(can_false && can_true);
  }
  if (answer)
    tree_add(x, y, true);
  else
    tree_add(y, x, false);
  return answer;
}

void adversary_decide(int i) { st.value[i] = st.decided++; }

// McIlroy-style lazy answers: every element starts undecided, and an undecided element is above
// every decided one. When two undecided elements meet, one of them is decided (the one last seen
// undecided, if it is one of the two), as the next value up. All answers are consistent with the
// values in the end. With tie_at >= 0, at that call every undecided element gets one of 1-4 values
// above all decided ones, so the rest of the sort runs on ties.
bool adversary_answer(int x, int y) {
  if (st.tie_at >= 0 && (long)st.answers.size() == st.tie_at) {
    int64_t top = st.decided;
    int64_t kinds = 1 + (int64_t)(st.tie_rng() % 4);
    for (int i = 0; i < st.n; i++)
      if (st.value[i] < 0) st.value[i] = top + (int64_t)(st.tie_rng() % kinds);
    st.decided = top + kinds;
  }
  if (st.value[x] < 0 && st.value[y] < 0) {
    if (x == st.candidate)
      adversary_decide(x);
    else
      adversary_decide(y);
  }
  if (st.value[x] < 0)
    st.candidate = x;
  else if (st.value[y] < 0)
    st.candidate = y;
  int64_t vx = st.value[x] < 0 ? INT64_MAX : st.value[x];
  int64_t vy = st.value[y] < 0 ? INT64_MAX : st.value[y];
  return vx < vy;
}

struct Compare {
  bool operator()(int x, int y) const {
    bool inside = x >= 0 && x < st.n && y >= 0 && y < st.n;
    bool answer = false;
    if (!inside) {
      st.oob++;
    } else {
      switch (st.mode) {
        case Mode::Keys: answer = st.key[x] < st.key[y]; break;
        case Mode::Ranges: answer = range_less(x, y); break;
        case Mode::Tree: answer = tree_answer(x, y); break;
        case Mode::Adversary: answer = adversary_answer(x, y); break;
      }
    }
    st.calls_a.push_back(x);
    st.calls_b.push_back(y);
    st.answers.push_back(answer);
    return answer;
  }
};

long total_oob = 0, total_sorts = 0, total_calls = 0;

// The elements sorted: plain ints (alg s, p), or a heavier type that carries the position with a
// string and a shared pointer, as NCBI's ranges with fuzz carry references (alg S, P). The
// comparator sees only the positions.
struct Heavy {
  int id;
  std::string pad;
  std::shared_ptr<int> ref;
};
int id_of(int e) { return e; }
int id_of(const Heavy& e) { return e.id; }
int make_elem(int id, int) { return id; }
Heavy make_elem(int id, Heavy) { return Heavy{id, std::string(24, 'a' + (id & 15)), std::make_shared<int>(id)}; }

template <typename E>
struct CompareE {
  bool operator()(const E& x, const E& y) const { return Compare{}(id_of(x), id_of(y)); }
};

// Sorts positions 0..n-1 with the current mode and writes the record line.
template <typename E>
void run_with(const char* tag, char alg, const std::string& extra) {
  int n = st.n;
  std::vector<E> buf(n + 2 * GUARD);
  for (int i = 0; i < GUARD; i++) {
    buf[i] = make_elem(-1 - i, E{});
    buf[GUARD + n + i] = make_elem(n + i, E{});
  }
  for (int i = 0; i < n; i++) buf[GUARD + i] = make_elem(i, E{});
  st.calls_a.clear();
  st.calls_b.clear();
  st.answers.clear();
  st.oob = 0;
  E* first = buf.data() + GUARD;
  E* last = first + n;
  if (alg == 's' || alg == 'S')
    std::sort(first, last, CompareE<E>{});
  else
    std::partial_sort(first, last, last, CompareE<E>{});
  for (int i = 0; i < GUARD; i++) {
    if (id_of(buf[i]) != -1 - i || id_of(buf[GUARD + n + i]) != n + i) st.oob++;  // a guard cell was written
  }
  total_oob += st.oob;
  total_sorts++;
  total_calls += (long)st.answers.size();
  std::string line;
  line.reserve(32 + st.answers.size() * 12 + (size_t)n * 6 + extra.size());
  line += tag;
  line += ' ';
  line += alg;
  line += ' ';
  line += std::to_string(n);
  line += ' ';
  if (st.answers.empty()) {
    line += "- -";
  } else {
    for (char a : st.answers) line += a ? '1' : '0';
    line += ' ';
    for (size_t k = 0; k < st.answers.size(); k++) {
      if (k) line += ',';
      line += std::to_string(st.calls_a[k]);
      line += ':';
      line += std::to_string(st.calls_b[k]);
    }
  }
  line += ' ';
  if (n == 0) {
    line += '-';
  } else {
    for (int i = 0; i < n; i++) {
      if (i) line += ',';
      line += std::to_string(id_of(first[i]));
    }
  }
  line += ' ';
  line += std::to_string(st.oob);
  if (!extra.empty()) {
    line += ' ';
    line += extra;
  }
  line += '\n';
  fwrite(line.data(), 1, line.size(), stdout);
}

void run_one(const char* tag, char alg, const std::string& extra) {
  if (alg == 'S' || alg == 'P')
    run_with<Heavy>(tag, alg, extra);
  else
    run_with<int>(tag, alg, extra);
}

std::string keys_text() {
  std::string s;
  for (int i = 0; i < st.n; i++) {
    if (i) s += ',';
    s += std::to_string(st.key[i]);
  }
  return s.empty() ? "-" : s;
}

void mode_tree(char alg, int n) {
  st.mode = Mode::Tree;
  st.n = n;
  st.tree_answers.clear();
  st.tree_alt.clear();
  st.tree_prefix = 0;
  long leaves = 0;
  for (;;) {
    st.le.assign((size_t)n * n, 0);
    st.lt.assign((size_t)n * n, 0);
    for (int i = 0; i < n; i++) st.le[(size_t)i * n + i] = 1;
    run_one("T", alg, "");
    leaves++;
    if (st.tree_answers.size() != st.answers.size()) {
      fprintf(stderr, "a run made fewer calls than its script\n");
      exit(2);
    }
    // The next leaf: the last call answered false whose answer true is consistent too.
    long k = (long)st.tree_answers.size() - 1;
    while (k >= 0 && !(st.tree_alt[k] && !st.tree_answers[k])) k--;
    if (k < 0) break;
    st.tree_answers.resize(k + 1);
    st.tree_alt.resize(k + 1);
    st.tree_answers[k] = 1;
    st.tree_alt[k] = 0;
    st.tree_prefix = k + 1;
  }
  fprintf(stderr, "tree %c n=%d leaves=%ld\n", alg, n, leaves);
}

// Lengths: uniform in 0..maxn, every fourth case 0..40, every eighth 15..70.
int pick_n(std::mt19937_64& g, long c, int maxn) {
  if (c % 8 == 3) return 15 + (int)(g() % 56);
  if (c % 4 == 1) return (int)(g() % 41);
  return (int)(g() % (uint64_t)(maxn + 1));
}

void mode_keys(char alg, uint64_t seed, long cases, int maxn) {
  st.mode = Mode::Keys;
  std::mt19937_64 g(seed);
  for (long c = 0; c < cases; c++) {
    int n = pick_n(g, c, maxn);
    st.n = n;
    st.key.assign(n, 0);
    int shape = (int)(g() % 12);
    uint64_t doms[] = {1, 2, 3, 4, 5, 8, 16, 64};
    uint64_t d;
    switch (shape) {
      case 0:  // small domain: many ties
      case 1:
        d = doms[g() % 8];
        for (auto& k : st.key) k = g() % d;
        break;
      case 2:  // domain about n/4 or n
        d = 1 + (uint64_t)n / (g() % 2 ? 4 : 1);
        for (auto& k : st.key) k = g() % d;
        break;
      case 3:  // distinct
        for (auto& k : st.key) k = g();
        break;
      case 4:  // ascending with ties
      case 5:  // descending with ties
        d = doms[g() % 8];
        for (auto& k : st.key) k = g() % d;
        std::sort(st.key.begin(), st.key.end());
        if (shape == 5) std::reverse(st.key.begin(), st.key.end());
        break;
      case 6:  // organ pipe
        for (int i = 0; i < n; i++) st.key[i] = (uint64_t)std::min(i, n - 1 - i) / (1 + g() % 3);
        break;
      case 7:  // saw
        d = 2 + g() % 20;
        for (int i = 0; i < n; i++) st.key[i] = (uint64_t)i % d;
        break;
      case 8:  // nearly sorted: a few swaps
        for (int i = 0; i < n; i++) st.key[i] = (uint64_t)i / (1 + g() % 2);
        for (int s = 0; n > 1 && s < 1 + (int)(g() % 4); s++)
          std::swap(st.key[g() % n], st.key[g() % n]);
        break;
      case 9:  // all equal
        break;
      case 10:  // two values alternating, or a block
        for (int i = 0; i < n; i++) st.key[i] = g() % 2 ? (uint64_t)(i % 2) : (uint64_t)(i < n / 2);
        break;
      default:  // mostly one value, a few others
        for (auto& k : st.key) k = g() % 10 == 0 ? g() % 50 : 7;
        break;
    }
    run_one("K", alg, std::to_string(shape) + "/" + keys_text());
  }
}

// NCBI-shaped ranges: empty ranges (to_open <= from, including an open end of 0), whole ranges,
// and normal ranges with duplicates and shared starts; the classes of the bounds set (exactly two
// different empty ranges, one empty range repeated, none) and random mixtures.
void mode_ranges(char alg, uint64_t seed, long cases, int maxn) {
  st.mode = Mode::Ranges;
  std::mt19937_64 g(seed);
  auto empty_range = [&](uint32_t& f, uint32_t& t) {
    switch (g() % 3) {
      case 0: f = (uint32_t)(g() % 10000); t = 0; break;  // to = 4294967295
      case 1: f = 0xFFFFFFFFu; t = 0xFFFFFFFFu; break;     // GetEmpty()
      default: f = 1 + (uint32_t)(g() % 10000); t = (uint32_t)(g() % f); break;
    }
  };
  auto normal = [&](uint32_t base, uint32_t& f, uint32_t& t) {
    if (g() % 5 == 0) {
      f = base + (uint32_t)(g() % 6);
      t = f + 1 + (uint32_t)(g() % 6);
    } else {
      f = (uint32_t)(g() % 10000);
      t = f + 1 + (uint32_t)(g() % 250);
    }
  };
  for (long c = 0; c < cases; c++) {
    int n = pick_n(g, c, maxn);
    st.n = n;
    st.reverse = g() % 2;
    st.from.assign(n, 0);
    st.to_open.assign(n, 0);
    uint32_t base = (uint32_t)(g() % 9000);
    int cls = (int)(g() % 5);
    for (int i = 0; i < n; i++) normal(base, st.from[i], st.to_open[i]);
    if (cls == 0 && n >= 2) {  // exactly two different empty ranges
      uint32_t f1, t1, f2, t2;
      empty_range(f1, t1);
      do empty_range(f2, t2); while (f1 == f2 && t1 == t2);
      int a = (int)(g() % n), b;
      do b = (int)(g() % n); while (b == a);
      st.from[a] = f1; st.to_open[a] = t1;
      st.from[b] = f2; st.to_open[b] = t2;
    } else if (cls == 1 && n >= 1) {  // one empty range repeated
      uint32_t f1, t1;
      empty_range(f1, t1);
      int k = 1 + (int)(g() % 4);
      for (int r = 0; r < k; r++) {
        int a = (int)(g() % n);
        st.from[a] = f1; st.to_open[a] = t1;
      }
    } else if (cls == 2) {  // none
    } else {  // a mixture: a share of empty ranges, some whole ones
      int pct = 1 + (int)(g() % 60);
      for (int i = 0; i < n; i++) {
        int x = (int)(g() % 100);
        if (x < pct) empty_range(st.from[i], st.to_open[i]);
        else if (x == pct && g() % 3 == 0) { st.from[i] = 0; st.to_open[i] = 0xFFFFFFFFu; }
      }
    }
    std::string extra = std::to_string(st.reverse ? 1 : 0) + "/";
    for (int i = 0; i < n; i++) {
      if (i) extra += ',';
      extra += std::to_string(st.from[i]) + ":" + std::to_string(st.to_open[i]);
    }
    if (n == 0) extra += '-';
    run_one("G", alg, extra);
  }
}

void mode_adv(char alg, uint64_t seed, long cases, int maxn) {
  st.mode = Mode::Adversary;
  std::mt19937_64 g(seed);
  for (long c = 0; c < cases; c++) {
    int n = pick_n(g, c, maxn);
    st.n = n;
    // First the plain adversary, to count its calls; then, every other case, a run that turns the
    // undecided elements into ties at a random call.
    for (int pass = 0; pass < 2; pass++) {
      st.value.assign(n, -1);
      st.decided = 0;
      st.candidate = -1;
      st.tie_at = -1;
      if (pass == 1) {
        if (c % 2 == 0) break;
        long calls = (long)st.answers.size();
        st.tie_at = calls ? (long)(g() % (uint64_t)calls) : 0;
        st.tie_rng.seed(g());
      }
      run_one("A", alg, std::to_string(st.tie_at));
    }
  }
}

void mode_file(char alg) {
  st.mode = Mode::Keys;
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    int n;
    if (!(in >> n)) continue;
    st.n = n;
    st.key.assign(n, 0);
    for (auto& k : st.key) in >> k;
    run_one("F", alg, "-/" + keys_text());
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: see the header of observe_sort.cpp\n");
    return 2;
  }
  std::string mode = argv[1];
  char alg = argv[2][0];
  if (alg != 's' && alg != 'p' && alg != 'S' && alg != 'P') {
    fprintf(stderr, "alg: s, p, S or P\n");
    return 2;
  }
  if (mode == "tree" && argc == 4) {
    mode_tree(alg, atoi(argv[3]));
  } else if ((mode == "keys" || mode == "ranges" || mode == "adv") && argc == 6) {
    uint64_t seed = strtoull(argv[3], nullptr, 0);
    long cases = atol(argv[4]);
    int maxn = atoi(argv[5]);
    if (mode == "keys") mode_keys(alg, seed, cases, maxn);
    else if (mode == "ranges") mode_ranges(alg, seed, cases, maxn);
    else mode_adv(alg, seed, cases, maxn);
  } else if (mode == "file") {
    mode_file(alg);
  } else {
    fprintf(stderr, "bad arguments\n");
    return 2;
  }
  fflush(stdout);
  fprintf(stderr, "sorts %ld calls %ld out-of-range %ld\n", total_sorts, total_calls, total_oob);
  return 0;
}
