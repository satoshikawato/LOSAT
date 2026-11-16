// The permutations of libstdc++'s std::sort (and std::partial_sort over the whole range) under
// NCBI's CRange_Less / CRange_ReverseLess on (from, to_open) ranges: six random cases, then each
// `CASE <reverse> <n> from:to_open ...` line read from the standard input (antiqsort.py).
#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
struct R { uint32_t f, t; };
static bool empty(R r){ return r.t <= r.f; }
static bool whole(R r){ return r.f==0 && r.t==0xFFFFFFFFu; }
static uint32_t to(R r){ return r.t-1; }
struct Less { bool rev; bool operator()(const R&a,const R&b) const {
  if (whole(a)) return !whole(b);
  if (empty(a)) return !empty(b) && !whole(b);
  if (whole(b)||empty(b)) return false;
  if (rev) { if (to(a)!=to(b)) return to(a)>to(b); return a.f<b.f; }
  return a.f<b.f || (a.f==b.f && a.t<b.t);
}};
static void print_case(bool rev, std::vector<R> v){
  printf("CASE %d %zu", rev, v.size());
  for(auto&r:v) printf(" %u:%u", r.f, r.t);
  std::vector<R> h=v;
  std::partial_sort(h.begin(), h.end(), h.end(), Less{rev});
  std::sort(v.begin(), v.end(), Less{rev});
  printf("\nSORTED");
  for(auto&r:v) printf(" %u:%u", r.f, r.t);
  printf("\nHEAP");
  for(auto&r:h) printf(" %u:%u", r.f, r.t);
  printf("\n");
}
int main(){
  std::mt19937 g(7);
  for (int c=0;c<6;c++){
    int n = c<3 ? 40+c*37 : 300+c;
    std::vector<R> v;
    for(int i=0;i<n;i++){
      int k=g()%5; R r;
      if(k==0){ r.f=g()%100+10; r.t=g()%r.f; }          // empty, distinct values
      else if(k==1){ r.f=g()%50; r.t=0; }                 // empty with open end 0
      else { r.f=g()%200; r.t=r.f+1+g()%30; }
      v.push_back(r);
    }
    print_case(c%2, v);
  }
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string word; int rev; size_t n;
    in >> word >> rev >> n;
    std::vector<R> v;
    for (size_t i = 0; i < n; i++) { R r; char colon; in >> r.f >> colon >> r.t; v.push_back(r); }
    print_case(rev, v);
  }
}
