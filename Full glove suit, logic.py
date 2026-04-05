import ui
SKADSU Developer API

Star Key Authentication Device Service Unit

This API treats SKADSU as one organism with five organs (HUB, ICU, COAT, Glove, RBP), exposed through a unified interface.

---

1. Concepts

• Star Key: Core identity + timing signature (lives on RBP).
• Tree: One of the 8 defense trees.
• Cadence Mode: DRIVE | FLOW | ANCHOR.
• States: IDLE, TREE_SELECTED, TREE_CONFIRMED, BLEND_MODE, UNIFIED_CADENCE, EIGHTFOLD_SHIELD.


---

2. Data Models

2.1 Tree

{
  "id": "ACTIVE",
  "finger": "L_INDEX"
}


Allowed: ACTIVE, PASSIVE, NEUTRAL, WATER, SPATIAL, EARTH, PLANETARY, TEAM.

2.2 CadenceMode

{
  "mode": "FLOW"
}


Allowed: DRIVE, FLOW, ANCHOR.

2.3 SKADSUState

{
  "state": "UNIFIED_CADENCE",
  "selectedTree": "ACTIVE",
  "secondaryTree": "SPATIAL",
  "cadenceMode": "FLOW",
  "unifiedCadenceActive": true,
  "eightfoldShieldActive": false
}


---

3. Transport

You can expose this over:

• BLE GATT,
• WebSocket, or
• local IPC.


Below is a logical API, not bound to a specific transport.

---

4. Commands

4.1 Get current SKADSU state

Command: GET_STATE

Request:

{ "cmd": "GET_STATE" }


Response:

{
  "cmd": "STATE",
  "payload": { /* SKADSUState */ }
}


---

4.2 Set cadence mode

Command: SET_CADENCE_MODE

Request:

{
  "cmd": "SET_CADENCE_MODE",
  "payload": { "mode": "ANCHOR" }
}


Effect:

• RBP adjusts timing curve
• HUB/ICU/COAT shift rhythm
• Glove palm haptics update


---

4.3 Select tree

Command: SELECT_TREE

Request:

{
  "cmd": "SELECT_TREE",
  "payload": { "tree": "WATER" }
}


Effect:

• Glove: TIP pattern on mapped finger
• RBP: load Star Key lane
• HUB/ICU/COAT: prep for that tree


---

4.4 Confirm tree

Command: CONFIRM_TREE

Request:

{
  "cmd": "CONFIRM_TREE"
}


Effect:

• Glove: MID pattern on selected finger
• RBP: generate timing signature
• COAT: open authority lane


---

4.5 Set blend

Command: SET_BLEND

Request:

{
  "cmd": "SET_BLEND",
  "payload": {
    "primary": "WATER",
    "secondary": "SPATIAL"
  }
}


Effect:

• Dual‐tree state
• ICU: dual HUD
• RBP: blended timing curve


---

4.6 Activate Unified Cadence

Command: ACTIVATE_UNIFIED_CADENCE

Request:

{
  "cmd": "ACTIVATE_UNIFIED_CADENCE"
}


Effect:

• RBP: Star Key resonance
• HUB: presence lock
• ICU: identity lock
• COAT: authority lock
• Glove: field haptics


---

4.7 Activate Eightfold Shield

Command: ACTIVATE_EIGHTFOLD_SHIELD

Request:

{
  "cmd": "ACTIVATE_EIGHTFOLD_SHIELD"
}


Precondition: all 8 trees seen this session.

Effect:

• RBP: sanctuary mode
• HUB: human‐first override
• ICU: full‐field visualization
• COAT: command sanctuary seal
• Glove: full cascade pattern


---

4.8 Quiet mode

Command: SET_QUIET_MODE

{
  "cmd": "SET_QUIET_MODE",
  "payload": { "enabled": true }
}


Effect: reduce haptic amplitude/duty cycle across all organs.

---

5. Events

Devices emit events upstream; apps can subscribe.

5.1 Tree changed

{
  "event": "TREE_CHANGED",
  "payload": {
    "selectedTree": "EARTH",
    "source": "GLOVE"
  }
}


5.2 Cadence mode changed

{
  "event": "CADENCE_MODE_CHANGED",
  "payload": {
    "mode": "DRIVE",
    "source": "HUB"
  }
}


5.3 Unified Cadence state

{
  "event": "UNIFIED_CADENCE_STATE",
  "payload": {
    "active": true,
    "source": "RBP"
  }
}


5.4 Eightfold Shield state

{
  "event": "EIGHTFOLD_SHIELD_STATE",
  "payload": {
    "active": false,
    "source": "RBP"
  }
}


---

6. Integration notes

• RBP is source of truth for identity & Star Key.
• Glove is primary input, but other organs may suggest mode changes.
• All devices should converge on the same SKADSUState within a small time window (e.g., <100ms on local network).
v = ui.load_view()
v.present('sheet')