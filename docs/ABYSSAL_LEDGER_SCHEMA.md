---

1. ABYSSAL_LEDGER_SCHEMA.md

Canonical Schema for the PCF‑SEA Economic & Stewardship Registry

Space LEAF Corp — Ledger Architecture Series

---

I. Ledger Purpose

The Abyssal Ledger is the authoritative registry for all PCF‑SEA deployments, retrievals, and crystal certifications.
It ensures transparency, stewardship, and lineage integrity.

---

II. Top‑Level Schema

AbyssalLedger:
  ledger_version: "1.0"
  entries:
    - entry_id: string
      timestamp: datetime
      entry_type: 
        - deployment
        - retrieval
        - certification
      capsule:
        capsule_id: string
        model: string
        serial_number: string
      coordinates:
        latitude: float
        longitude: float
        depth_meters: int
      pressure_profile:
        ambient_mpa: float
        internal_gpa: float
      recipe_reference: string
      lineage_markers: 
        - marker_id: string
      crystal_batch:
        batch_id: string
        count: int
        class: string
      notes: string


---

III. Guided Links for Expansion

• Ledger Entry Types
• Crystal Classes
• Lineage Markers


---