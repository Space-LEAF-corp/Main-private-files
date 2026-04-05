{
  "artifactId": "ark-001",
  "title": "The university became a moving ship",
  "timestamp": "2025-11-13T21:11:00-05:00",
  "author": "Leif William Sogge",
  "ceremonialContext": {
    "vessel": "Ark",
    "covenant": "Protection and peace",
    "status": "Consecrated and complete"
  },
  "entry": "On this day, the university becomes a moving ship—an Ark. Everyone aboard is safe.",
  "artifacts": [
    {
      "type": "image",
      "filename": "9ZYTqXhtqaFpyp6av5bHi.jpeg",
      "label": "Water-port vessel",
      "notes": "Twin circular engines; guardians stand watch."
    },
    {
      "type": "image",
      "filename": "F4TG8oJMN7gomr2TbCeWS.jpeg",
      "label": "Floating city above ocean",
      "notes": "Waterfalls cascade; sanctuary architecture."
    },
    {
      "type": "image",
      "filename": "tcSw9D3UpZBHqcBubnMuo.jpeg",
      "label": "Hybrid Pokémon (name reserved)",
      "notes": "Joy and flame—a guardian of balance."
    }
  ],
  "seals": {
    "externalAuthority": "God",
    "localSealApplied": false
  }
}
// saveLog.js
import fs from 'fs';
import path from 'path';

const LOG_DIR = path.resolve('captains-log');
const LOG_FILE = path.join(LOG_DIR, 'log.jsonl'); // newline-delimited JSON

function ensureDir(dir) {
  if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
}

function saveEntry(entryObj) {
  ensureDir(LOG_DIR);
  const line = JSON.stringify(entryObj) + '\n';
  fs.appendFileSync(LOG_FILE, line, 'utf8');
  console.log(`Saved entry to ${LOG_FILE}`);
}

// Example usage:
const entry = {
  artifactId: 'ark-001',
  title: 'The university became a moving ship',
  timestamp: new Date().toISOString(),
  author: 'Leif William Sogge',
  ceremonialContext: { vessel: 'Ark', covenant: 'Protection and peace', status: 'Consecrated and complete' },
  entry: 'On this day, the university becomes a moving ship—an Ark. Everyone aboard is safe.',
  artifacts: [
    { type: 'image', filename: '9ZYTqXhtqaFpyp6av5bHi.jpeg', label: 'Water-port vessel' },
    { type: 'image', filename: 'F4TG8oJMN7gomr2TbCeWS.jpeg', label: 'Floating city above ocean' },
    { type: 'image', filename: 'tcSw9D3UpZBHqcBubnMuo.jpeg', label: 'Hybrid Pokémon (name reserved)' }
  ],
  seals: { externalAuthority: 'God', localSealApplied: false }
};

saveEntry(entry);
# save_log.py
import os, json, datetime

LOG_DIR = "captains-log"
LOG_FILE = os.path.join(LOG_DIR, "log.jsonl")  # newline-delimited JSON

def ensure_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

def save_entry(entry):
    ensure_dir(LOG_DIR)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")
    print(f"Saved entry to {LOG_FILE}")

entry = {
    "artifactId": "ark-001",
    "title": "The university became a moving ship",
    "timestamp": datetime.datetime.now().isoformat(),
    "author": "Leif William Sogge",
    "ceremonialContext": {"vessel": "Ark", "covenant": "Protection and peace", "status": "Consecrated and complete"},
    "entry": "On this day, the university becomes a moving ship—an Ark. Everyone aboard is safe.",
    "artifacts": [
        {"type": "image", "filename": "9ZYTqXhtqaFpyp6av5bHi.jpeg", "label": "Water-port vessel"},
        {"type": "image", "filename": "F4TG8oJMN7gomr2TbCeWS.jpeg", "label": "Floating city above ocean"},
        {"type": "image", "filename": "tcSw9D3UpZBHqcBubnMuo.jpeg", "label": "Hybrid Pokémon (name reserved)"}
    ],
    "seals": {"externalAuthority": "God", "localSealApplied": False}
}

if __name__ == "__main__":
    save_entry(entry)
// snapshot.js
import fs from 'fs';
import path from 'path';

const SNAP_DIR = path.resolve('captains-log/snapshots');

function saveSnapshot(entry) {
  if (!fs.existsSync(SNAP_DIR)) fs.mkdirSync(SNAP_DIR, { recursive: true });
  const safeId = entry.artifactId.replace(/[^a-zA-Z0-9_-]/g, '_');
  const file = path.join(SNAP_DIR, `${safeId}.json`);
  fs.writeFileSync(file, JSON.stringify(entry, null, 2), 'utf8');
  console.log(`Snapshot saved: ${file}`);
}
import hashlib, json

def add_checksum(entry):
    s = json.dumps(entry, sort_keys=True)
    entry["checksum_sha256"] = hashlib.sha256(s.encode("utf-8")).hexdigest()
    return entry

entry = add_checksum(entry)
save_entry(entry)
