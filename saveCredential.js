// saveCredential.js
// Node.js 18+
// Purpose: Create + sign a consent-gated QR-DNA manifest and write to an append-only, hash-chained ledger.

const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

// ---- Configuration (set paths and identities) ----
const LEDGER_DIR = path.join(__dirname, 'ledger');
const LEDGER_FILE = path.join(LEDGER_DIR, 'qr_dna_ledger.jsonl');

// Master Steward private key (PEM, Ed25519 recommended). Store offline in production.
const MASTER_PRIV_KEY_PEM = process.env.MASTER_PRIV_KEY_PEM; // set via env var
const COUNCIL_PRIV_KEYS_PEM = (process.env.COUNCIL_PRIV_KEYS_PEM || '')
  .split('||')
  .filter(Boolean); // e.g., "-----BEGIN PRIVATE KEY-----...||-----BEGIN PRIVATE KEY-----..."

// Utility: load previous ledger head (hash) for chain linkage.
function getPrevHead() {
  if (!fs.existsSync(LEDGER_FILE)) return null;
  const fd = fs.openSync(LEDGER_FILE, 'r');
  try {
    // Read last line efficiently
    const stat = fs.fstatSync(fd);
    const buffer = Buffer.alloc(Math.min(stat.size, 1024 * 128));
    const start = Math.max(0, stat.size - buffer.length);
    fs.readSync(fd, buffer, 0, buffer.length, start);
    const lines = buffer.toString('utf8').trim().split('\n').filter(Boolean);
    const last = lines[lines.length - 1];
    if (!last) return null;
    const entry = JSON.parse(last);
    return entry.head; // previous chain head
  } catch {
    return null;
  } finally {
    fs.closeSync(fd);
  }
}

// Utility: base64url
function b64url(buf) {
  return Buffer.from(buf)
    .toString('base64')
    .replace(/\+/g, '-')
    .replace(/\//g, '_')
    .replace(/=+$/g, '');
}

// ---- Consent token (ephemeral, subject-signed in production) ----
function createConsentToken({ subjectId, scope, expiresAt }) {
  // Minimal unsigned consent token; in production, have the subject sign this with their device key.
  return {
    ver: 1,
    subjectId,           // pseudonymous ID, not raw PII
    scope,               // e.g., ["age_proof", "membership_check"]
    issuedAt: new Date().toISOString(),
    expiresAt,           // ISO time, short window
    nonce: b64url(crypto.randomBytes(24)),
  };
}

// ---- Create credential manifest ----
function createManifest({ stewardId, leaderId, jurisdiction, claims, consentToken }) {
  return {
    schema: 'ark.qr-dna.manifest.v1',
    stewardId,                 // your public steward ID
    leaderId,                  // head-of-state or delegate public ID
    jurisdiction,              // ISO 3166 or diplomatic code
    claims,                    // requested zero-knowledge claims (no raw DNA)
    consentToken,              // ephemeral consent artifact
    createdAt: new Date().toISOString(),
    artifactId: b64url(crypto.randomBytes(32)), // locally unique ID
  };
}

// ---- Signing helpers (Ed25519) ----
function signObject(obj, privPem, signerId) {
  const data = Buffer.from(JSON.stringify(obj));
  const key = crypto.createPrivateKey(privPem);
  const sig = crypto.sign(null, data, key); // Ed25519: pass null for algorithm with key type
  return {
    signerId, // public identifier for who signed
    alg: 'Ed25519',
    sig: b64url(sig),
  };
}

function computeHead(entry) {
  // Hash the entry deterministically for chain head
  const h = crypto.createHash('sha256');
  h.update(Buffer.from(JSON.stringify(entry)));
  return b64url(h.digest());
}

// ---- Append-only ledger write ----
function appendToLedger(entry) {
  if (!fs.existsSync(LEDGER_DIR)) fs.mkdirSync(LEDGER_DIR, { recursive: true });
  const line = JSON.stringify(entry) + '\n';
  fs.appendFileSync(LEDGER_FILE, line, { encoding: 'utf8', mode: 0o600 });
}

// ---- Main save function ----
function saveCredential({
  stewardId,
  leaderId,
  jurisdiction,
  claims,
  subjectId,
  consentScope,
  consentTTLMinutes = 10,
}) {
  if (!MASTER_PRIV_KEY_PEM) throw new Error('MASTER_PRIV_KEY_PEM not set');
  if (COUNCIL_PRIV_KEYS_PEM.length === 0) throw new Error('COUNCIL_PRIV_KEYS_PEM empty');

  // 1) Consent token
  const expiresAt = new Date(Date.now() + consentTTLMinutes * 60_000).toISOString();
  const consentToken = createConsentToken({ subjectId, scope: consentScope, expiresAt });

  // 2) Manifest (no personal payloads; only claims)
  const manifest = createManifest({
    stewardId,
    leaderId,
    jurisdiction,
    claims,
    consentToken,
  });

  // 3) Signatures (Master + Council quorum)
  const masterSig = signObject(manifest, MASTER_PRIV_KEY_PEM, stewardId);
  const councilSigs = COUNCIL_PRIV_KEYS_PEM.map((pem, idx) =>
    signObject(manifest, pem, `council-${idx + 1}`)
  );

  // 4) Build ledger entry (minimal disclosure + chain)
  const prevHead = getPrevHead();
  const entry = {
    type: 'qr-dna.verification.request.v1',
    manifest,
    signatures: {
      master: masterSig,
      council: councilSigs, // threshold policy is enforced by readers/verifiers
    },
    prevHead: prevHead || null,
  };

  // Compute head and append
  entry.head = computeHead(entry);
  appendToLedger(entry);

  // 5) Return minimal receipt (no payloads)
  return {
    receiptId: manifest.artifactId,
    head: entry.head,
    prevHead: entry.prevHead,
    createdAt: manifest.createdAt,
    jurisdiction,
    // Only pass/fail should be shared after verification; here we just confirm save.
    status: 'SAVED',
  };
}

// ---- Example usage (run with env vars set) ----
if (require.main === module) {
  const receipt = saveCredential({
    stewardId: 'ark-steward:leif',
    leaderId: 'state-head:example-country',
    jurisdiction: 'US',
    claims: ['age_proof', 'membership_check'], // zero-knowledge claims to be proven later
    subjectId: 'subject:anon-12345',
    consentScope: ['age_proof'],
    consentTTLMinutes: 5,
  });
  console.log('Saved:', receipt);
}
