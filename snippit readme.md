DEMO EASTER EGG — NOT FOR PRODUCTION

Purpose: playful, auditable demo to prove authorship and transparency.
Scope: staging/local only. This override is intentionally isolated and cannot access production secrets.
Activation: requires a signed short-lived token and explicit developer approval.
Kill switch: set ENV DEMO_OVERRIDE=false or run ./scripts/kill_override.sh
Audit: all attempts logged to /var/log/demo_override.audit (append-only, checksummed).
