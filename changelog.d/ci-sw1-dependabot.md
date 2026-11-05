- Dependencies: Dependabot version updates (`.github/dependabot.yml`) open weekly pull requests
  against `dev` for cargo (`LOSAT/`; before merging, `web/adapter/Cargo.lock` is moved to the same
  versions, as `web/adapter/tools/check_build_identity.py` requires), npm (`web/app/`; minor and patch
  updates in one pull request) and GitHub Actions, at most three open per ecosystem, for releases at
  least seven days old; Auto-merge merges them as it does every other pull request.
  `.github/CODEOWNERS` names the owner of seqdb, writedb, `.github/`, AGENTS.md and the product
  decisions.
