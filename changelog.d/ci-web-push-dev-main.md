- CI: web.yml runs on pushes only to main and dev, and by hand (`gh workflow run web.yml --ref
  <branch>`). A push to another branch ran the whole Web set before its pull request existed, and
  the pull request's run tested it again (42 % of the repository's job-minutes on 2026-11-11..12).
