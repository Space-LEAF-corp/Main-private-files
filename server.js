// server.js — private minimal backup for chalk bundles
// Run: node server.js
import express from "express";
import fs from "fs";
import path from "path";

const app = express();
app.use(express.json({ limit: "5mb" }));

const DATA_DIR = path.resolve("chalk_backups");
if (!fs.existsSync(DATA_DIR)) fs.mkdirSync(DATA_DIR);

app.post("/backup", (req, res) => {
  const bundle = req.body;
  if (bundle?.type !== "chalk-bundle" || !bundle?.integrity) {
    return res.status(400).json({ ok: false, error: "Invalid bundle" });
  }
  try {
    const stamp = new Date().toISOString().replace(/[:.]/g,"-");
    const file = path.join(DATA_DIR, `eternal-chalkboard-${stamp}.chalk.json`);
    fs.writeFileSync(file, JSON.stringify(bundle, null, 2), "utf8");
    return res.json({ ok: true, file });
  } catch (error) {
    return res.status(500).json({ ok: false, error: "Failed to write backup file" });
  }
});

app.listen(8787, () => console.log("Backup server on http://localhost:8787"));
