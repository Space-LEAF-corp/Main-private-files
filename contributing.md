### CONTRIBUTING.md

# Contributing to Scroll of Universal Respect 1.0

Thank you for your interest in contributing. This repository prioritizes consent, boundary safety, and auditability. Contributions must preserve those principles. The project is currently staged privately for partner review and security validation; public contributions will be welcomed after partner sign‑offs.

---

### Purpose

This document explains how to contribute, the review and release flow, code standards, testing expectations, and how security reports are handled by our partner hosts: Space Leaf Corp, Microsoft, and Tesla (silent monitoring partner).

---

### Getting Started

**Access and setup**
- Request access from the maintainers if the repository is private.
- Clone the repository and create a feature branch using the pattern `feat/<short-description>` or `fix/<short-description>`.
- Install dependencies and run the test suite locally before opening a pull request.

**Branching and commits**
- Branch prefixes: `feat/`, `fix/`, `chore/`, `docs/`, `test/`.
- Commit message style: imperative and scoped. Example:  
  `feat(guard): add consent token TTL enforcement`
- Signed commits: GPG or SSH-signed commits are required for protected branches.

---

### Code Standards and Tests

**Languages and style**
- Primary language: Python 3.11+.
- Follow PEP 8 for Python. Use `ruff` or `flake8` for linting.

**Testing**
- All new features must include unit tests.
- Tests must run with `pytest`.
- CI runs tests on every PR and push. Ensure your branch passes CI before requesting review.

**Static analysis and security**
- Run static analysis and dependency checks locally (e.g., `ruff`, `bandit`, `safety`).
- Never commit secrets, keys, or credentials. Use environment variables or GitHub Secrets for CI.

---

### Pull Request Process

**Before opening a PR**
- Rebase or merge the latest `main` into your branch.
- Run `pytest -q` and linters locally.
- Update or add documentation for new behavior.
- Add or update tests to cover edge cases and failure modes.

**PR content**
- Title: concise and prefixed with type (feat/fix/docs/test).
- Description: include motivation, what changed, and any migration steps.
- Include this checklist in the PR body:
  - [ ] Tests added or updated
  - [ ] Lint passed locally
  - [ ] Docs updated if behavior changed
  - [ ] Security considerations reviewed
  - [ ] Commits signed

**Review and merging**
- PRs require at least one approving review and passing CI checks.
- Address requested changes promptly.
- Maintain a clean history: squash or rebase as requested by maintainers.
- Merges to `main` are protected and require status checks and signed commits.

---

### Security and Responsible Disclosure

**Private staging**
- The repository is currently private for partner review. Do not publish or disclose private details externally.

**Reporting vulnerabilities**
- Do not open public issues for security matters.
- Follow the instructions in `SECURITY.md` to report vulnerabilities privately to the joint intake addresses.
- Microsoft and Space Leaf Corp co-host intake and remediation. Tesla participates as a silent monitoring partner for agreed telemetry and incident correlation under confidentiality agreements.

**Security expectations for contributors**
- Avoid introducing features that simulate intimacy or bypass consent checks.
- Add audit logging for any change that affects consent, policy rules, or telemetry.
- Document telemetry or monitoring changes and obtain legal/security sign-off before enabling sharing with partners.

---

### Governance Release and Roles

**Release gating**
- Public release requires:
  - Passing CI and security scans
  - Legal and MOU sign-offs for telemetry sharing with Tesla
  - Replacement of placeholder PGP fingerprints and contact emails in `SECURITY.md`
  - Presence of `LICENSE`, `CODE_OF_CONDUCT.md`, and `CONTRIBUTING.md`

**Roles**
- Maintainers: review PRs, merge after checks and approvals, manage releases.
- Security leads: triage advisories, coordinate with Microsoft and Tesla per `SECURITY.md`.
- Legal contacts: finalize MOUs and data processing agreements before telemetry sharing.

**Release process**
- Tag releases using semantic versioning (for example, `v1.0.0`).
- Include release notes summarizing changes, security considerations, and partner acknowledgements.
- Ensure signed release artifacts where applicable.

---

### Templates and Quick Snippets

**PR checklist**
