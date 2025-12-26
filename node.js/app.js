require("dotenv").config();
const express = require("express");
const bodyParser = require("body-parser");
const {
  hashPassword,
  checkPassword,
  generateToken,
  authMiddleware
} = require("./auth");

const app = express();
app.use(bodyParser.json());

// In-memory store for demo; replace with a real DB.
const users = [];

// Bootstrap a single admin on first run.
async function bootstrapAdmin() {
  const existing = users.find(u => u.role === "admin");
  if (existing) return;

  const email = process.env.ADMIN_EMAIL;
  const password = process.env.ADMIN_INITIAL_PASSWORD;

  if (!email || !password) {
    console.error("ADMIN_EMAIL and ADMIN_INITIAL_PASSWORD must be set");
    process.exit(1);
  }

  const passwordHash = await hashPassword(password);

  users.push({
    id: 1,
    email,
    passwordHash,
    role: "admin"
  });

  console.log("Admin user bootstrapped with configured credentials.");
  console.log("Remember to change ADMIN_INITIAL_PASSWORD and rotate it regularly.");
}

app.post("/auth/login", async (req, res) => {
  const { email, password } = req.body;
  const user = users.find(u => u.email === email);
  if (!user) {
    return res.status(401).json({ error: "Invalid credentials" });
  }

  const ok = await checkPassword(password, user.passwordHash);
  if (!ok) {
    return res.status(401).json({ error: "Invalid credentials" });
  }

  // MFA hooks would go here (TOTP, push, etc.)
  const token = generateToken(user);
  res.json({ token });
});

// Example: regular user route
app.get("/me", authMiddleware(), (req, res) => {
  res.json({ user: req.user });
});

// Example: admin-only sensitive route
app.get("/admin/status", authMiddleware("admin"), (req, res) => {
  res.json({
    status: "ok",
    message: "Admin route reached. Geometry intact.",
    timestamp: new Date().toISOString()
  });
});

const port = process.env.APP_PORT || 8080;
bootstrapAdmin().then(() => {
  app.listen(port, () => {
    console.log(`Space Leaf Platform listening on port ${port}`);
  });
});
