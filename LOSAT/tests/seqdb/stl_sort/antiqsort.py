#!/usr/bin/env python3
"""Search an input that drives libstdc++'s introsort (std::sort) into its depth-limit heap fallback
with empty ranges in the part it heap-sorts. Model of __introsort_loop / __unguarded_partition_pivot
(as LOSAT's range_set::stl_sort), with McIlroy's adversary ("A Killer Adversary for Quicksort",
1999) mirrored: undecided ("gas") elements compare below decided ones, and an element is decided
with the next value from the top. The adversary stops at the heap fallback: the elements still
undecided there are only known to be below every decided one, so they may be any such ranges,
empty ones among them. Prints `CASE <reverse> <n> from:to_open ...` lines for sort_cases.cpp.
"""
import sys

THRESHOLD = 16


class Adversary:
    def __init__(self, n):
        self.val = [None] * n  # None = gas
        self.top = n
        self.candidate = 0
        self.stopped = False

    def freeze(self, x):
        self.top -= 1
        self.val[x] = self.top

    def less(self, x, y):
        v = self.val
        if not self.stopped:
            if v[x] is None and v[y] is None:
                self.freeze(x if x == self.candidate else y)
            if v[x] is None:
                self.candidate = x
            elif v[y] is None:
                self.candidate = y
        a = -1 if v[x] is None else v[x]
        b = -1 if v[y] is None else v[y]
        return a < b


def introsort(v, less, on_fallback):
    n = len(v)
    depth = 2 * (n.bit_length() - 1)

    def loop(first, last, depth):
        while last - first > THRESHOLD:
            if depth == 0:
                on_fallback(first, last)
                return True
            depth -= 1
            mid = first + (last - first) // 2
            a, b, c = first + 1, mid, last - 1
            if less(v[a], v[b]):
                pick = b if less(v[b], v[c]) else (c if less(v[a], v[c]) else a)
            elif less(v[a], v[c]):
                pick = a
            elif less(v[b], v[c]):
                pick = c
            else:
                pick = b
            v[first], v[pick] = v[pick], v[first]
            lo, hi = first + 1, last
            while True:
                while less(v[lo], v[first]):
                    lo += 1
                hi -= 1
                while less(v[first], v[hi]):
                    hi -= 1
                if lo >= hi:
                    break
                v[lo], v[hi] = v[hi], v[lo]
                lo += 1
            if loop(lo, last, depth):
                return True
            last = lo
        return False

    return loop(0, n, depth)


def make(n, reverse, empties_every):
    adv = Adversary(n)
    order = list(range(n))
    seen = {}

    def on_fallback(first, last):
        adv.stopped = True
        seen["range"] = (first, last, list(order[first:last]))

    if not introsort(order, adv.less, on_fallback):
        return None
    first, last, part = seen["range"]
    gas = [i for i in range(n) if adv.val[i] is None]
    # Decided elements: non-empty ranges above every undecided one, in the order of their values.
    ranges = [None] * n
    for i in range(n):
        v = adv.val[i]
        if v is None:
            continue
        if reverse:
            # CRange_ReverseLess: a larger `to` comes first; values grow as `to` falls.
            ranges[i] = (v, 1000 - v + 1)  # to = 1000 - v
        else:
            ranges[i] = (1000 + v, 1000 + v + 3)
    # Undecided elements: empty ranges (all equivalent, different) and small non-empty ones.
    for k, i in enumerate(gas):
        if k % empties_every == 0:
            ranges[i] = (40 + k % 37, k % 29)  # empty: to_open <= from
        elif reverse:
            ranges[i] = (k % 50, 2000 + (k * 7) % 300)  # to above every decided `to`
        else:
            ranges[i] = (k * 7 % 500, k * 7 % 500 + 1 + k % 13)
    return ranges, (first, last), len(gas), sum(1 for k in range(len(gas)) if k % empties_every == 0)


def main():
    for reverse, n in ((0, 120), (1, 150)):
        made = make(n, reverse, 3)
        if made is None:
            sys.exit(f"no fallback for n={n}")
        ranges, (first, last), gas, empties = made
        print(f"n={n} reverse={reverse} fallback on [{first},{last}) gas={gas} empties={empties}",
              file=sys.stderr)
        print(f"CASE {reverse} {n} " + " ".join(f"{a}:{b}" for a, b in ranges))


if __name__ == "__main__":
    main()
