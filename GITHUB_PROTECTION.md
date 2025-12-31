GitHub Branch Protection & Required Reviewers
===========================================

These are suggested steps an admin can run or perform via the GitHub UI to enforce the governance policies declared in `GOVERNANCE.md`.

Via the GitHub UI (recommended for non-automated setup):

1. Go to `Settings -> Branches` for the repository.
2. Under "Branch protection rules", add a rule for branch `main`.
3. Enable "Require pull request reviews before merging" and set the required number of approving reviews (e.g., 1 or 2).
4. Enable "Require review from Code Owners" so that entries in `.github/CODEOWNERS` are enforced.
5. Optionally enable "Require status checks to pass before merging" and add CI checks.

Via `gh` (GitHub CLI) — example (requires admin privileges):

```bash
# install gh and authenticate as an admin
gh auth login

# Create branch protection requiring code owner review (example)
gh api --method PUT /repos/:owner/:repo/branches/main/protection -f required_status_checks=null -f enforce_admins=true -f required_pull_request_reviews='{"dismiss_stale_reviews":true,"require_code_owner_reviews":true,"required_approving_review_count":1}'
```

Notes

- The `require_code_owner_reviews` flag makes GitHub require review from the users/teams listed in `.github/CODEOWNERS` for affected files.
- Only a repository admin can apply these settings.
