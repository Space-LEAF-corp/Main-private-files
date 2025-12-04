GOVERNANCE
==========

Purpose
-------
This document records the governance policy for core laws and settings in this repository.

Authority
---------
- Leif William Sogge (referred to as "Leif") is designated as the sole code owner authorized to alter "core laws" and core settings without explicit UN approval.
- Changes to core laws or settings must be made via a pull request and require review by Leif. For additional security, Microsoft review is also requested via CODEOWNERS.

UN Oversight
-----------
- Where the repository license or external agreements require United Nations (UN) approval for certain changes, UN approval must be obtained prior to merging such changes. This document does not supersede binding external agreements; it only records repository-level policy.

Procedure for core changes
--------------------------
1. Open a pull request targeting `main` and reference the related governance discussion or UN approval record when applicable.
2. Assign `@GuardianNinja` (Leif) as a required reviewer. The PR must include a clear justification and link to any UN approval if applicable.
3. Microsoft (`@microsoft`) is also requested as a reviewer via `.github/CODEOWNERS` and should give a security / partner review.
4. After required approvals, an admin may merge the PR. If UN approval is required by external governance, a record of such approval must be included in the PR description.

Contact
-------
- Security contact: `leifwsogge@me.com` (per `SECURITY.md`)

Notes
-----
- This is a repository-level governance file. Enforcing these rules for branches, required reviewers and status checks requires repository admin settings (branch protection rules). See `GITHUB_PROTECTION.md` for suggested admin steps.
