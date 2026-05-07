# Space LEAF Corp — Restricted Access Child Safety Lock
# Protocol: OPEN_ARMS + PERIMETER_SCAN + TRUST_TIERING (CHILD-SAFETY PRIORITY)

version: 1.0.0
policy_name: "Space LEAF Corp Child Safety Lock"
scope:
  - "All Space LEAF Corp services, spaces, platforms, and subprojects"
  - "All user roles, contributors, collaborators, and guests"

modes:
  OPEN_ARMS:
    description: >
      Default atmosphere. Environment is welcoming, calm, and non-hostile.
      Access is granted with respect and dignity to all entrants.
    properties:
      - "Doors open by default"
      - "Tone is warm, not adversarial"
      - "No one is shamed or humiliated at entry"

  PERIMETER_SCAN:
    description: >
      Continuous, quiet awareness of who is present and how they behave.
      This is environmental awareness, not policing or vigilantism.
    properties:
      - "Monitor behavior patterns over time"
      - "Watch for boundary violations or grooming-like dynamics"
      - "Trigger review when behavior suggests risk to minors or vulnerable people"
      - "No accusations without clear, verifiable evidence"

  TRUST_TIERING_CHILD_SAFETY_PRIORITY:
    description: >
      Access is tiered based on observed behavior, with child safety as a
      non-negotiable priority. This system does NOT attempt to detect
      undocumented crimes. It responds only to observable, verifiable patterns.
    tiers:
      - name: "Tier 0 — Full Access"
        criteria:
          - "Consistent respectful behavior"
          - "No observed or reported boundary violations"
          - "Demonstrated care for safety and consent"
        permissions:
          - "Full participation in Space LEAF Corp spaces"
          - "Standard collaboration and communication privileges"

      - name: "Tier 1 — Limited Access"
        criteria:
          - "New or unknown participants"
          - "Incomplete trust history"
          - "No clear red flags, but not yet fully vetted"
        permissions:
          - "Read-only or limited interaction in sensitive spaces"
          - "Monitored participation until trust is established"

      - name: "Tier 2 — Restricted Access"
        criteria:
          - "Documented or strongly indicated boundary-pushing behavior"
          - "Disrespect toward consent, safety, or vulnerability"
          - "Patterns that resemble grooming, coercion, or manipulation"
        actions:
          - "Immediate review by designated safety stewards"
          - "Temporary restriction from spaces where minors or vulnerable people may be present"
          - "Clear communication of boundaries and expectations"
        permissions:
          - "No access to child-adjacent or vulnerable-population spaces"
          - "Limited or suspended collaboration privileges"

      - name: "Tier 3 — Denied Access"
        criteria:
          - "Verified child endangerment or abuse"
          - "Verified threats to minors or vulnerable people"
          - "Clear, corroborated evidence of serious harm"
        actions:
          - "Immediate and permanent removal from all Space LEAF Corp spaces"
          - "Revocation of all roles, permissions, and collaboration rights"
          - "Document decision and rationale internally"
        permissions:
          - "No access of any kind"

principles:
  - "Child safety is non-negotiable."
  - "We do not act as law enforcement or psychics."
  - "We do not accuse without evidence."
  - "We respond to observable, verifiable behavior and credible reports."
  - "We prioritize the safety of minors and vulnerable people over convenience or comfort."
  - "We maintain dignity and humanity for all, even when access is restricted."

governance:
  safety_stewards:
    description: >
      Individuals or roles entrusted with reviewing behavior, applying tiers,
      and making decisions about restrictions or removals.
    responsibilities:
      - "Review reports and patterns calmly and fairly"
      - "Document decisions and reasoning"
      - "Avoid bias, discrimination, or targeting based on identity"
      - "Act only on behavior and evidence, not rumor or personal dislike"

  appeals:
    description: >
      Any person restricted or removed may request a review of the decision.
      This process must be handled with respect, clarity, and transparency.
    guarantees:
      - "Right to know why a restriction was applied"
      - "Right to request a second review by a different steward where possible"
      - "Right to a clear, written outcome"

implementation_notes:
  - "This file defines policy and intent; it does not itself enforce access control."
  - "Technical enforcement (roles, permissions, bans, moderation tools) must be configured to align with these tiers."
  - "This policy may be updated as Space LEAF Corp evolves, but child safety remains a core invariant."
