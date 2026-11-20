#!/usr/bin/env bash
# Is the commit of this run the source of the release's tag?
#
# usage: release_tag_source.sh <version> <true|false: does this run upload>
#
# A publishing run builds the commit it was started on; the release's tag is checked here and never
# checked out (LOSAT/tests/test_prepare_release_candidate_v020.py: no workflow checks out an input
# ref). The check passes when this commit is the tag v<version>, or differs from it only under
# .github/ (a fix of the publishing workflows after the tag publishes the same source). Anything
# else is an error in a run that uploads, and a notice in a rehearsal.
set -euo pipefail

version=$1
uploads=$2
tag="v$version"

not_the_tag() {
  if [ "$uploads" = true ]; then
    echo "::error::$1"
    exit 1
  fi
  echo "::notice::rehearsal: $1"
  exit 0
}

if ! git fetch --quiet --no-tags --depth=1 origin "+refs/tags/$tag:refs/tags/$tag"; then
  not_the_tag "there is no tag $tag, and this commit's version is $version"
fi
if ! git diff --quiet "refs/tags/$tag" HEAD -- . ':(exclude).github'; then
  not_the_tag "this commit ($(git rev-parse HEAD)) is not the source of the tag $tag: start the run on the tag, or on a branch that differs from it only under .github/"
fi
echo "The source of this run is the tag $tag ($(git rev-parse "refs/tags/$tag^{commit}"))."
