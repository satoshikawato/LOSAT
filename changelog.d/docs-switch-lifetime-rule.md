- AGENTS.md states the life of a speed-up switch (`LOSAT_X_*`): it lands on by default in its own pull request,
  marks its rule with `x_reach!`, has a unit test against its `=0` function and an output test against `=0`, has
  no `SHADOW` / `VERIFY` mode, and graduates at the first graduation check after it reaches `main`.
