const bcrypt = require("bcrypt");
const jwt = require("jsonwebtoken");

const SALT_ROUNDS = 12; // can increase if CPU allows

async function hashPassword(password) {
  return bcrypt.hash(password, SALT_ROUNDS);
}

async function checkPassword(password, hash) {
  return bcrypt.compare(password, hash);
}

function generateToken(user) {
  // Never put secrets like passwords in the token.
  const payload = {
    id: user.id,
    role: user.role, // "admin" or "user"
    email: user.email
  };
  return jwt.sign(payload, process.env.JWT_SECRET, {
    expiresIn: "1h"
  });
}

function authMiddleware(requiredRole) {
  return function (req, res, next) {
    const authHeader = req.headers.authorization || "";
    const token = authHeader.replace("Bearer ", "").trim();
    if (!token) {
      return res.status(401).json({ error: "Missing token" });
    }

    try {
      const decoded = jwt.verify(token, process.env.JWT_SECRET);
      req.user = decoded;

      if (requiredRole && decoded.role !== requiredRole) {
        return res.status(403).json({ error: "Forbidden" });
      }

      next();
    } catch (err) {
      return res.status(401).json({ error: "Invalid token" });
    }
  };
}

module.exports = {
  hashPassword,
  checkPassword,
  generateToken,
  authMiddleware
};
