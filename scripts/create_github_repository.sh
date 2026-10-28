#!/usr/bin/env bash
set -euo pipefail

# Use the user's existing gh authentication. Never read or write tokens here.
expected_owner=FacundoFirmenich
expected_name=crispr-observability
expected_repository="$expected_owner/$expected_name"
repository_root=$(git rev-parse --show-toplevel)
cd "$repository_root"
command -v gh >/dev/null || { echo 'GitHub CLI (gh) is required.' >&2; exit 1; }
active_owner=$(gh api user --jq .login)
[[ "$active_owner" == "$expected_owner" ]] || { echo 'Active GitHub account differs from FacundoFirmenich.' >&2; exit 1; }
git rev-parse --verify HEAD >/dev/null
[[ -z "$(git status --porcelain)" ]] || { echo 'Commit or remove working-tree changes before publication.' >&2; exit 1; }
[[ "$(git branch --show-current)" == main ]] || { echo 'Expected local branch main.' >&2; exit 1; }
[[ -z "$(git remote)" ]] || { echo 'Existing Git remotes require explicit inspection; no changes made.' >&2; exit 1; }
if gh repo view "$expected_repository" --json nameWithOwner >/dev/null 2>&1; then
  echo 'Target repository already exists; inspect it before any upload. No changes made.' >&2
  exit 1
fi
# An unavailable repository lookup can also mean an API/network failure. The
# create operation below is authoritative and refuses an occupied repository.
gh repo create "$expected_repository" --public --source=. --remote=origin --push \
  --description 'Product-dependent observability and sequence correspondence in genome editing readouts: public figures, numerical data and computational checks'
remote_url=$(git remote get-url origin)
case "$remote_url" in
  "https://github.com/$expected_repository.git"|"git@github.com:$expected_repository.git"|"https://github.com/$expected_repository") ;;
  *) echo 'Unexpected resulting origin URL; inspect the repository before continuing.' >&2; exit 1 ;;
esac
gh repo view "$expected_repository" --json nameWithOwner,url,visibility,defaultBranchRef
