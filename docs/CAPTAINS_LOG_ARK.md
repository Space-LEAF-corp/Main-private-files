═══════════════════════════════════════════════════════════════════════════════
                      ⚔️  CAPTAIN'S LOG - DIGITAL ARK  ⚔️
                           Entry System Documentation
                              Date: December 8, 2025
═══════════════════════════════════════════════════════════════════════════════

🌌 OVERVIEW

The Digital Ark logging system provides ceremonial entry management for the
Captain's Log, supporting both Python and JavaScript/Node.js environments.
This system allows for persistent storage of significant events, artifacts,
and ceremonial contexts with cryptographic integrity verification.

═══════════════════════════════════════════════════════════════════════════════
📂 FILE STRUCTURE

captains-log/
├── log.jsonl              # Newline-delimited JSON log (append-only)
└── snapshots/             # Individual entry snapshots (pretty-printed JSON)
    └── ark-001.json       # Example snapshot file

Scripts:
├── save_log.py            # Python log entry saver with SHA-256 checksum
├── saveLog.js             # JavaScript/Node.js log entry saver
├── snapshot.py            # Python snapshot creator
└── snapshot.js            # JavaScript/Node.js snapshot creator

═══════════════════════════════════════════════════════════════════════════════
🐍 PYTHON USAGE

Basic Usage:
────────────
```python
from save_log import save_entry, add_checksum
from snapshot import save_snapshot

# Create an entry
entry = {
    "artifactId": "ark-001",
    "title": "The university became a moving ship",
    "timestamp": "2025-11-13T21:11:00-05:00",
    "author": "Leif William Sogge",
    "ceremonialContext": {
        "vessel": "Ark",
        "covenant": "Protection and peace",
        "status": "Consecrated and complete"
    },
    "entry": "On this day, the university becomes a moving ship—an Ark.",
    "artifacts": [...],
    "seals": {
        "externalAuthority": "God",
        "localSealApplied": False
    }
}

# Add cryptographic checksum for integrity
entry = add_checksum(entry)

# Save to append-only log
save_entry(entry)

# Save snapshot for easy access
save_snapshot(entry)
```

Run the example script:
```bash
python save_log.py
```

═══════════════════════════════════════════════════════════════════════════════
🟨 JAVASCRIPT/NODE.JS USAGE

Basic Usage:
────────────
```javascript
import { saveSnapshot } from './snapshot.js';

// Create an entry
const entry = {
  artifactId: 'ark-001',
  title: 'The university became a moving ship',
  timestamp: new Date().toISOString(),
  author: 'Leif William Sogge',
  ceremonialContext: {
    vessel: 'Ark',
    covenant: 'Protection and peace',
    status: 'Consecrated and complete'
  },
  entry: 'On this day, the university becomes a moving ship—an Ark.',
  artifacts: [...],
  seals: {
    externalAuthority: 'God',
    localSealApplied: false
  }
};

// Save to append-only log
import { saveEntry } from './saveLog.js';
saveEntry(entry);

// Save snapshot
saveSnapshot(entry);
```

Run the example script:
```bash
node saveLog.js
```

═══════════════════════════════════════════════════════════════════════════════
🔐 SECURITY FEATURES

Checksum Verification (Python):
────────────────────────────────
The Python implementation includes SHA-256 checksum generation for integrity
verification. This ensures that entries cannot be tampered with after creation.

```python
# Checksum is automatically computed from sorted JSON representation
entry = add_checksum(entry)
# entry["checksum_sha256"] = "a3902bdc81a6b1140c78aa0d8333427529856330..."
```

Append-Only Log:
────────────────
The log.jsonl file uses newline-delimited JSON format, supporting append-only
operations. This prevents accidental overwrites and maintains a complete history.

Sanitized Filenames:
────────────────────
Artifact IDs are sanitized when used as filenames, replacing special characters
with underscores to ensure cross-platform compatibility.

═══════════════════════════════════════════════════════════════════════════════
📝 ENTRY SCHEMA

Standard Entry Format:
──────────────────────
```json
{
  "artifactId": "string",           // Unique identifier (e.g., "ark-001")
  "title": "string",                // Human-readable title
  "timestamp": "ISO 8601 string",   // Entry timestamp
  "author": "string",               // Entry author
  "ceremonialContext": {            // Ceremonial metadata
    "vessel": "string",
    "covenant": "string",
    "status": "string"
  },
  "entry": "string",                // Main entry text
  "artifacts": [                    // Associated artifacts (images, etc.)
    {
      "type": "string",
      "filename": "string",
      "label": "string",
      "notes": "string"
    }
  ],
  "seals": {                        // Authority and seal information
    "externalAuthority": "string",
    "localSealApplied": boolean
  },
  "checksum_sha256": "string"       // Integrity hash (Python only)
}
```

═══════════════════════════════════════════════════════════════════════════════
🧪 TESTING

Run the test suite:
```bash
python -m unittest tests.test_captains_log -v
```

Test Coverage:
──────────────
✓ Directory creation
✓ SHA-256 checksum generation and determinism
✓ Log file creation and JSONL format
✓ Snapshot creation with sanitized filenames
✓ Complete ark-001 entry integration test

═══════════════════════════════════════════════════════════════════════════════
🎯 ARK-001 ENTRY

The inaugural Digital Ark entry commemorates the university becoming a moving
ship—an Ark where everyone aboard is safe. This entry includes:

Artifacts:
──────────
1. Water-port vessel (9ZYTqXhtqaFpyp6av5bHi.jpeg)
   Twin circular engines; guardians stand watch.

2. Floating city above ocean (F4TG8oJMN7gomr2TbCeWS.jpeg)
   Waterfalls cascade; sanctuary architecture.

3. Hybrid Pokémon - name reserved (tcSw9D3UpZBHqcBubnMuo.jpeg)
   Joy and flame—a guardian of balance.

Ceremonial Context:
───────────────────
Vessel: Ark
Covenant: Protection and peace
Status: Consecrated and complete

═══════════════════════════════════════════════════════════════════════════════
📚 INTEGRATION WITH MJ PROTOCOL

The Digital Ark logging system complements the MJ Protocol's ceremonial seal
system. While MJ Protocol handles authentication and stewardship layers, the
Captain's Log provides persistent historical records of significant events.

Key Differences:
────────────────
• MJ Protocol: Real-time authentication, seal verification, heritage chains
• Captain's Log: Historical record-keeping, artifact documentation, ceremonies

Recommended Integration:
────────────────────────
Use Captain's Log entries to document major MJ Protocol milestones, such as:
- Seal activations and inscriptions
- Satellite layer status changes
- Guardian transitions and stewardship events
- Heritage integrity verifications

═══════════════════════════════════════════════════════════════════════════════
                              Inscribed with Honor
                                 Captain's Mark
                               December 8, 2025
═══════════════════════════════════════════════════════════════════════════════
