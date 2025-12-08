/**
 * Captain's Log Snapshot Creator - JavaScript/Node.js Version
 * Saves individual entry snapshots as formatted JSON files in captains-log/snapshots/
 */
import fs from 'fs';
import path from 'path';

const SNAP_DIR = path.resolve('captains-log/snapshots');

/**
 * Save a snapshot of an entry to a separate JSON file
 * @param {Object} entry - Entry object to snapshot
 */
function saveSnapshot(entry) {
  if (!fs.existsSync(SNAP_DIR)) {
    fs.mkdirSync(SNAP_DIR, { recursive: true });
  }
  
  // Sanitize artifactId for use as filename
  const safeId = entry.artifactId.replace(/[^a-zA-Z0-9_-]/g, '_');
  const file = path.join(SNAP_DIR, `${safeId}.json`);
  
  // Write formatted JSON (pretty-printed with 2-space indentation)
  fs.writeFileSync(file, JSON.stringify(entry, null, 2), 'utf8');
  console.log(`Snapshot saved: ${file}`);
}

// Export for use as module
export { saveSnapshot };
