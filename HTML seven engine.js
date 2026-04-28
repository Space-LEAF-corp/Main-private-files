
---

## 🧩 How Modes Work

Modes are JSON files that define:

- emoji sets  
- animation style + speed  
- audio behavior  
- background  
- density  
- visual effects  

See `/docs/mode-schema.md` for the full specification.

---

## 🛠️ Running the Project

1. Clone the repo  
2. Open `index.html` in a browser  
3. Switch modes by editing:

```js
import { MODES } from "./modes.js";
engine.loadMode(MODES.Calm);
