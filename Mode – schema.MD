
---

# 📘 **mode-schema.md (the “HTML Seven” language spec)**

```markdown
# Mode Schema Specification (HTML Seven)

Each mode is a JSON object describing how the engine should render emoji
textures, animations, and audio behavior.

---

## 🔧 Schema

```json
{
  "emojiSet": ["✨", "🔥"],
  "animation": {
    "style": "pulse",
    "speed": 1.0
  },
  "audioBehavior": {
    "trigger": "click",
    "autoplay": false
  },
  "background": "linear-gradient(#000, #111)",
  "density": 100,
  "effects": {
    "glow": false,
    "drift": true
  }
}
