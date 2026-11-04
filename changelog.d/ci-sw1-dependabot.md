- Dependencies: Dependabot version updates (`.github/dependabot.yml`) open weekly pull requests
  against `dev` for cargo (LOSAT and web/adapter together, one pull request per dependency, so that
  both lock files keep the same versions), npm (web/app; minor and patch updates in one pull
  request) and GitHub Actions, at most three open per ecosystem, for releases at least seven days
  old. Auto-merge leaves these pull requests to a person. `.github/CODEOWNERS` names the owner of
  seqdb, writedb, `.github/`, AGENTS.md and the product decisions.
