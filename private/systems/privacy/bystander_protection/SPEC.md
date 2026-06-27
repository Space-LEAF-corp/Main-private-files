---

BYSTANDER PROTECTION MODULE (BPM)

Space LEAF Corp — Internal Specification

File: /private/systems/privacy/bystander_protection/SPEC.md
Version: 1.0
Classification: INTERNAL — DO NOT DISTRIBUTE

---

1. Overview

The Bystander Protection Module (BPM) ensures that no non‑consenting individual becomes identifiable “content” in any media captured by Space LEAF Corp smart glasses or related devices. All detection, anonymization, and enforcement occur locally, with strict context‑aware rules and no raw unprotected frames ever written to storage.

This module integrates with:

• Consent Ledger
• Guardian Upload Gatekeeper
• Kid‑Safe Auto‑Blur Engine


---

2. Core Principles

2.1 Local‑First Processing

All detection, classification, and anonymization occur on‑device.
No cloud inference.
No external transmission of raw frames.

2.2 Default Anonymization

Any face not explicitly tagged as a primary subject is treated as a bystander and anonymized.

2.3 Context‑Aware Enforcement

Stricter rules apply in sensitive environments:

• Schools
• Hospitals
• Gyms
• Homes
• Worship spaces
• Courtrooms
• Changing areas


---

3. System Pipeline

3.1 Capture Layer

• Input: Raw video frames from smart glasses.
• Action: Frames are immediately routed to the Bystander Detector before any storage or buffering.
• Guarantee: No unprocessed frame is ever written to disk.


---

3.2 Bystander Detector

3.2.1 Face Detection

Function:

f(\text{frame}) \rightarrow \{\text{face\_boxes}\}


3.2.2 Primary Subject Classification

A face is considered primary if:

• It matches the user’s own face (self‑capture).
• It is explicitly tagged as a consenting subject (via Consent Ledger).


3.2.3 Bystander Rule

Any face not classified as primary → BYSTANDER.

---

3.3 Protection Actions

3.3.1 Real‑Time Blur

All bystander regions are anonymized using:

• Gaussian blur
• Pixelation
• Or full masking (context‑dependent)


3.3.2 No Raw Frame Storage

Only the protected frame is ever written to disk.

3.3.3 Metadata Flag

Frames containing bystanders must include:

BYSTANDER_PRESENT = true


---

4. User Controls

4.1 Guardian Mode (Default)

• All bystanders blurred.
• No override permitted.
• Mandatory in all sensitive contexts.


4.2 Ceremonial Mode

• Allows temporary unblur only if a matching Consent Ledger entry exists.
• Cannot override child protection rules.


4.3 HUD Indicator

A simple icon must display:
“Bystanders Detected — Protected.”

---

5. Context‑Aware Rules

5.1 Sensitive Contexts

The BPM must enforce stricter anonymization in:

• Hospitals / ICU
• Schools
• Gyms
• Homes
• Worship spaces
• Courtrooms
• Changing rooms
• Pools
• Locker rooms


5.2 Hard‑Block Zones

In the following contexts, recording is allowed but uploading is permanently blocked, regardless of consent:

• Hospitals / ICU
• Changing rooms
• Locker rooms
• Therapy offices
• Confessionals
• Private counseling spaces


---

6. Integration Points

6.1 Consent Ledger

The BPM queries the ledger to determine:

• Whether a face is a primary subject
• Whether unblur is permitted in Ceremonial Mode


6.2 Guardian Upload Gatekeeper

The BPM provides:

• BYSTANDER_PRESENT flag
• Context tags
• Protected frame metadata


The Gatekeeper uses this to decide:

• Allow upload
• Block upload
• Require consent
• Force Private Vault storage


6.3 Kid‑Safe Auto‑Blur Engine

If a face is classified as CHILD, the BPM defers to the stricter child‑protection rules.

---

7. Performance Requirements

• Latency: Must operate in real‑time on wearable hardware.
• Accuracy:• Face detection ≥ 95% in normal lighting
• Misclassification of children must err on the side of safety

• Power Efficiency:• Optimized for low‑power inference
• No unnecessary cloud calls



---

8. Security Requirements

• All processing must occur locally.
• No biometric templates stored.
• No raw frames cached.
• All metadata must be signed by the device key.
• No external service may request unprotected frames.


---

9. Testing Requirements

9.1 Unit Tests

• Face detection accuracy
• Primary subject classification
• Blur application
• Metadata flagging


9.2 Scenario Tests

• Crowded public spaces
• Schools
• Hospitals
• Gyms
• Homes
• Low‑light environments
• Rapid motion


9.3 Regression Tests

• Ensure no raw frame storage
• Ensure blur cannot be bypassed
• Ensure context rules remain enforced


---

10. Compliance & Audit

• All BPM actions must be logged in the Guardian Audit Log.
• Logs must contain no biometric or identifying data.
• Logs must be immutable and signed.


---

END OF SPECIFICATION

Space LEAF Corp — Internal Use Only

---
