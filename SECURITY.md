# SECURITY.md

## Purpose
This document explains how to report security vulnerabilities and how Microsoft and Space Leaf Corp operate jointly as Prime Cyber Security Services to receive, triage, and remediate reports for this repository.

## Scope
This policy covers:
- repository code and configuration
- published packages and releases
- CI CD pipelines and GitHub Actions
- infrastructure as code and deployment manifests

## Reporting a Vulnerability
To report a security issue, follow these steps:
1. Create a private GitHub Security Advisory for this repository with a clear title and reproduction steps.
2. Email the joint intake inboxes **security@spaceleafcorp.com** and **security@microsoft.com** with the subject line: Security Report: Scroll of Universal Respect.
3. If you prefer encrypted reports, use the PGP keys listed below.
4. Do not publish details publicly until the issue is resolved and disclosure is coordinated.

## PGP Keys
- **Space Leaf Corp PGP fingerprint** 0123 4567 89AB CDEF 0123 4567 89AB CDEF 0123 45 (placeholder — replace with real key)
- **Microsoft security PGP fingerprint** 89AB CDEF 0123 4567 89AB CDEF 0123 4567 89AB CD (placeholder — replace with official key)

## What to Include in a Report
Please include:
- a concise summary of the issue and potential impact
- step by step reproduction steps and a minimal test case
- affected versions, commit SHAs, and environment details
- suggested mitigations or patches if available
- contact information for follow up and preferred disclosure name or anonymity request

## Response Process and Timelines
- **Acknowledgement within 48 hours** of receipt
- **Initial triage and severity assessment within 5 business days**
- **Mitigation plan or patch timeline communicated within 10 business days**
- **Coordinated disclosure and public advisory** after fixes are verified

## Confidentiality and Safe Harbor
- Reports will be treated confidentially until coordinated disclosure
- Researchers acting in good faith who follow this policy will not be pursued for legal action, subject to applicable law and the researcher’s adherence to this policy

## Disclosure and Credits
- Reporters will be credited by name or handle unless they request anonymity
- Public advisories will include impact, root cause, and remediation steps

## Contact and Escalation
- **Primary intake** security@spaceleafcorp.com
- **Co host intake** security@microsoft.com
- If no response within 72 hours, escalate to the repository owners and organization security contacts

## Additional Notes
- Do not include sensitive data such as private keys or passwords in public issues or pull requests
- Use private advisories or encrypted email for sensitive attachments
- Replace placeholder PGP fingerprints with real keys before publishing this file
