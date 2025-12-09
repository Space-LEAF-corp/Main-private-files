// server.js — private minimal backup for chalk bundles
// Run: node server.js
// 
// SECURITY NOTE: This is a minimal backup server for private use.
// For production deployment, consider adding:
// - Rate limiting (e.g., express-rate-limit) to prevent DoS attacks
// - Authentication/authorization for the /backup endpoint
// - Disk space monitoring to prevent filling up storage
// - HTTPS/TLS for encrypted communication
import express from "express";
import fs from "fs";
import path from "path";

const app = express();
app.use(express.json({ limit: "5mb" }));

const DATA_DIR = path.resolve("chalk_backups");
fs.mkdirSync(DATA_DIR, { recursive: true });

app.post("/backup", (req, res) => {
  const bundle = req.body;
  // Validate bundle is a plain object, not null, not array
  if (
    typeof bundle !== "object" ||
    bundle === null ||
    Array.isArray(bundle)
  ) {
    return res.status(400).json({ ok: false, error: "Bundle must be a plain object" });
  }
  // Validate required fields and types
  if (
    typeof bundle.type !== "string" ||
    bundle.type !== "chalk-bundle" ||
    typeof bundle.integrity !== "string"
  ) {
    return res.status(400).json({ ok: false, error: "Invalid bundle fields" });
  }
  // Only allow expected properties (type, integrity, data)
  const allowedProps = ["type", "integrity", "data"];
  for (const key of Object.keys(bundle)) {
    if (!allowedProps.includes(key)) {
      return res.status(400).json({ ok: false, error: "Unexpected bundle property: " + key });
    }
  }
  // Check serialized size (max 5mb, but be strict)
  const jsonStr = JSON.stringify(bundle, null, 2);
  if (Buffer.byteLength(jsonStr, "utf8") > 5 * 1024 * 1024) {
    return res.status(413).json({ ok: false, error: "Bundle too large" });
  }
  try {
    const stamp = new Date().toISOString().replace(/[:.]/g, "-");
    const file = path.join(DATA_DIR, `eternal-chalkboard-${stamp}.chalk.json`);
    fs.writeFileSync(file, jsonStr, "utf8");
    return res.json({ ok: true, file });
  } catch (error) {
    console.error("Failed to write backup file:", error);
    return res.status(500).json({ ok: false, error: "Failed to write backup file" });
  }
});

app.listen(8787, () => console.log("Backup server on http://localhost:8787"));
