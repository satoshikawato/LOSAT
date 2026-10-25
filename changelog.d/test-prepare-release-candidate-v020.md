- Release process: `LOSAT/tests/prepare_release_candidate_v020.py` assembles the eight v0.2.0
  release-candidate artifacts (five native, serial and threaded Wasm, the source crate) and the
  handoff, after checking the candidate against the contract: HEAD is the candidate, the certified
  commit `C` is an ancestor, no runtime input differs between `C` and the candidate, the certification
  runs succeeded, and the evidence hashes are filled. Each artifact is smoked once per program with
  the NCBI-certified hashes. The v0.1.0 tooling is unchanged. No change to `losat`.
