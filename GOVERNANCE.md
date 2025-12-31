GOVERNANCE
==========

This document records the governance policy for core laws and settings in this repository.

Purpose
-------

This document records the governance policy for core laws and settings in this repository.

Authority
---------

- Leif William Sogge (referred to as "Leif") is designated as a primary code owner for core laws and core settings.
- Changes to core laws or settings require the explicit approval of Leif plus dual approval from both the President of the United States (placeholder GitHub handle `@POTUS-placeholder`) and `@GuardianNinja` (the repository owner). In other words, all three approvals are required to merge changes classified as "core." Replace `@POTUS-placeholder` with the official GitHub account or team representing the President of the United States or an authorized representative.
- All changes must be made via a pull request and include links to any external approvals when applicable.

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
