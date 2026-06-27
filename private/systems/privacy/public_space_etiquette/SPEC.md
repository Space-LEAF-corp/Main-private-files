---

PUBLIC‑SPACE ETIQUETTE PROTOCOL (PSEP)

Space LEAF Corp — Internal Specification

File: /private/systems/privacy/public_space_etiquette/SPEC.md
Version: 1.0
Classification: INTERNAL — DO NOT DISTRIBUTE

---

1. Overview

The Public‑Space Etiquette Protocol (PSEP) defines the behavioral, visual, and social rules that govern how Space LEAF Corp smart glasses operate in public environments. Its purpose is to ensure that the presence of the device never creates fear, confusion, or unwanted exposure, and that all interactions remain respectful, transparent, and dignity‑preserving.

PSEP is the human‑facing layer of the Privacy Suite.
It governs how the device signals recording, how it responds to bystander concerns, and how users are guided to behave ethically.

---

2. Core Principles

2.1 Transparency

People must always know when recording is active.

2.2 Respect for Personal Boundaries

If someone expresses discomfort, the system must respond immediately and respectfully.

2.3 No Stealth Recording

The device must never record secretly in public spaces.

2.4 Context Sensitivity

Certain environments require stricter etiquette and stronger protections.

2.5 Human‑Centered Design

The protocol must be intuitive, predictable, and socially considerate.

---

3. Recording Signals

3.1 Visual Indicator

A visible LED or icon on the glasses must illuminate when recording is active.

Requirements:

• Must be clearly visible from the front.
• Must not be user‑disableable.
• Must activate instantly when capture begins.


---

3.2 Audio Indicator (Optional)

A soft chime may play when recording starts or stops.

Requirements:

• Must be subtle, non‑alarming.
• Must not be used in sensitive environments (hospitals, therapy offices).


---

3.3 HUD Status

The user must always see a clear HUD indicator showing:

• LIVE CAPTURE: ON
• LIVE CAPTURE: OFF
• BYSTANDERS PROTECTED


---

4. Respect Zones

4.1 Hard No‑Record Zones

Recording is allowed, but upload is permanently blocked, regardless of consent:

• Bathrooms
• Changing rooms
• Locker rooms
• Therapy offices
• Confessionals
• Private counseling spaces


Requirements:

• Device must display:
“Sensitive Zone — Upload Disabled.”


---

4.2 Soft Consent Zones

Recording allowed only with explicit consent or event‑level permission:

• Gyms
• Schools
• Workplaces
• Worship spaces


Requirements:

• Device must prompt user:
“Consent required in this environment.”


---

4.3 Public Spaces

Recording allowed with standard protections:

• Parks
• Streets
• Public transit
• Outdoor events


Requirements:

• Bystander Protection Module active
• Child blur active
• Recording signal visible


---

5. Approach Protocol

5.1 If someone asks, “Are you recording?”

The device must display a clear, immediate status screen:

• LIVE CAPTURE: OFF
• or
• LIVE CAPTURE: ON (BYSTANDERS PROTECTED)


Requirements:

• Must be accessible within 1 second.
• Must not reveal private user data.


---

5.2 If someone says, “Please don’t record me.”

The system must:

1. Mark the individual as Do‑Not‑Capture for the session.
2. Apply auto‑blur to all future appearances.
3. Confirm visually to the user:
“Do‑Not‑Capture active.”


Requirements:

• No override permitted.
• Applies even in Ceremonial Mode.


---

5.3 If someone expresses discomfort non‑verbally

Examples:

• Covering face
• Turning away
• Raising hand


The system must treat this as a Do‑Not‑Capture signal.

---

6. Social Scripts (User Guidance)

The device must provide optional, user‑friendly scripts to help users communicate respectfully.

Approved Scripts:

• “These glasses blur everyone around me by default.”
• “I’m only recording myself; you’re protected.”
• “If you’d like, I can show you how it anonymizes you.”
• “I can turn off recording if you prefer.”


Requirements:

• Scripts must be accessible via quick HUD menu.
• Scripts must never imply that privacy is optional.


---

7. Public‑Facing Charter

The PSEP must include a public document titled:

“Space LEAF Corp Smart Glasses Etiquette Charter”

This document must explain:

• No stealth recording
• Bystander protection
• Kid‑safe rules
• Sensitive zone restrictions
• Upload limitations
• Respect‑first philosophy


Requirements:

• Must be written in plain language
• Must be accessible on the public website
• Must not reveal internal implementation details


---

8. Integration Points

8.1 Bystander Protection Module

PSEP relies on BPM for:

• Bystander detection
• Auto‑blur
• Do‑Not‑Capture enforcement


---

8.2 Kid‑Safe Auto‑Blur Engine

PSEP enforces:

• Mandatory child protection
• No override in sensitive zones


---

8.3 Guardian Upload Gatekeeper

PSEP defines:

• Context rules
• Sensitive zone restrictions
• Consent requirements


Gatekeeper enforces upload decisions.

---

9. Security Requirements

• Recording signals must not be spoofable.
• Do‑Not‑Capture flags must be cryptographically enforced.
• Sensitive zone detection must be tamper‑resistant.
• No user setting may disable etiquette protections.


---

10. Testing Requirements

10.1 Unit Tests

• LED activation
• HUD status accuracy
• Do‑Not‑Capture flagging
• Sensitive zone detection


10.2 Scenario Tests

• Public park
• Gym
• School
• Hospital
• Worship space
• Workplace
• Crowded environments


10.3 Regression Tests

• Ensure recording signal cannot be disabled
• Ensure Do‑Not‑Capture cannot be bypassed
• Ensure sensitive zone rules remain enforced


---

END OF SPECIFICATION

Space LEAF Corp — Internal Use Only

---
