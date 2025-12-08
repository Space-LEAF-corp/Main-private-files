/**
 * Captain's Log Entry Saver - JavaScript/Node.js Version
 * Saves ceremonial entries to captains-log/log.jsonl in newline-delimited JSON format.
 */
import fs from 'fs';
import path from 'path';

const LOG_DIR = path.resolve('captains-log');
const LOG_FILE = path.join(LOG_DIR, 'log.jsonl'); // newline-delimited JSON

/**
 * Ensure directory exists, create if not.
 * @param {string} dir - Directory path
 */
function ensureDir(dir) {
  if (!fs.existsSync(dir)) {
    fs.mkdirSync(dir, { recursive: true });
  }
}

/**
 * Save entry to log.jsonl
 * @param {Object} entryObj - Entry object to save
 */
function saveEntry(entryObj) {
  ensureDir(LOG_DIR);
  const line = JSON.stringify(entryObj) + '\n';
  fs.appendFileSync(LOG_FILE, line, 'utf8');
  console.log(`Saved entry to ${LOG_FILE}`);
}

// Export for use as module
export { saveEntry };

// Example usage when run directly: ark-001 - The university became a moving ship
if (import.meta.url === `file://${process.argv[1]}`) {
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
  entry: 'On this day, the university becomes a moving ship—an Ark. Everyone aboard is safe.',
  artifacts: [
    {
      type: 'image',
      filename: '9ZYTqXhtqaFpyp6av5bHi.jpeg',
      label: 'Water-port vessel',
      notes: 'Twin circular engines; guardians stand watch.'
    },
    {
      type: 'image',
      filename: 'F4TG8oJMN7gomr2TbCeWS.jpeg',
      label: 'Floating city above ocean',
      notes: 'Waterfalls cascade; sanctuary architecture.'
    },
    {
      type: 'image',
      filename: 'tcSw9D3UpZBHqcBubnMuo.jpeg',
      label: 'Hybrid Pokémon (name reserved)',
      notes: 'Joy and flame—a guardian of balance.'
    }
  ],
  seals: {
    externalAuthority: 'God',
    localSealApplied: false
  }
  };

  saveEntry(entry);
}
