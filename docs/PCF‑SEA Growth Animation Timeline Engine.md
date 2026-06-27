---

🌊 PCF‑SEA Growth Animation Timeline Engine

Abyssal Time‑Lapse of Crystal Thickening

---

I. Growth Model (Pressure‑Driven Expansion)

Crystal radius \( R(t) \) grows according to:

R(t) = R_0 + \alpha \cdot P_{\text{eff}} \cdot S \cdot t


Where:

• \( R_0 \) = initial seed radius
• \( \alpha \) = growth coefficient derived from K
• \( P_{\text{eff}} \) = effective pressure (GPa)
• \( S \) = supersaturation (0–1)
• \( t \) = time in months


We compute:

\alpha = \frac{1}{K}


So:

• High K → slow growth
• Low K → fast growth


This ties directly into your existing recipe system.

---

II. Timeline Simulation (Months → Frames)

We simulate 12 months by default:

• 1 month = 30 frames
• 12 months = 360 frames


Each frame:

• Expands lattice spacing
• Adds new atoms
• Extends bonds outward
• Updates 3D projection


---

III. Growth Animation Engine (JS Module)

Add this to your webapp’s <script> section:

let growthFrame = 0;
let growthTimeline = [];
let growthAnimating = false;

function generateGrowthTimeline(months = 12) {
  const K = parseFloat(document.getElementById("K_out").innerText);
  const P_eff = parseFloat(document.getElementById("P_eff_out").innerText);
  const S = parseFloat(document.getElementById("S_val").innerText);

  const alpha = 1 / K;
  const frames = months * 30;

  growthTimeline = [];

  for (let f = 0; f < frames; f++) {
    const t = f / 30; // months
    const scale = 1 + alpha * P_eff * S * t;
    growthTimeline.push(scale);
  }
}

function startGrowthAnimation() {
  growthAnimating = true;
  growthFrame = 0;
  animateGrowth();
}

function stopGrowthAnimation() {
  growthAnimating = false;
}

function animateGrowth() {
  if (!growthAnimating) return;

  const scale = growthTimeline[growthFrame] || 1;
  drawMoleculeGrowth(scale);

  growthFrame++;
  if (growthFrame < growthTimeline.length) {
    requestAnimationFrame(animateGrowth);
  } else {
    growthAnimating = false;
  }
}

function drawMoleculeGrowth(scale) {
  const canvas = document.getElementById("moleculeCanvas");
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  const cx = canvas.width / 2;
  const cy = canvas.height / 2;

  const baseScale = 60 * scale;

  const projected = moleculeData.atoms.map(a => {
    const x = cx + a.x * baseScale;
    const y = cy - a.y * baseScale - a.z * 10;
    return { x, y, element: a.element };
  });

  ctx.strokeStyle = "#555";
  ctx.lineWidth = 2;
  moleculeData.bonds.forEach(b => {
    const a1 = projected[b.a];
    const a2 = projected[b.b];
    ctx.beginPath();
    ctx.moveTo(a1.x, a1.y);
    ctx.lineTo(a2.x, a2.y);
    ctx.stroke();
  });

  projected.forEach(p => {
    ctx.beginPath();
    ctx.fillStyle = atomColor(p.element);
    ctx.arc(p.x, p.y, 8, 0, 2 * Math.PI);
    ctx.fill();
  });
}


---

IV. UI Controls

Add these buttons under your 3D viewer:

<button onclick="generateGrowthTimeline(12)">Generate 12‑Month Timeline</button>
<button onclick="startGrowthAnimation()">Play Growth Animation</button>
<button onclick="stopGrowthAnimation()">Stop</button>


---

V. How It Looks in Motion

The animation shows:

• Seed molecule →
• Lattice expansion →
• Bonds stretching →
• New layers forming →
• Crystal thickening over time


It becomes a living pressure‑driven organism, growing under the ocean’s weight.

---

VI. Integration With Recipes

When a recipe is selected:

loadSeedForClass(recipe.suggestedClass);
generateGrowthTimeline(recipe.estimatedMonths);
startGrowthAnimation();


This means:

• Every recipe has a seed molecule
• Every seed has a growth timeline
• Every timeline becomes a visual story of its formation


This is the ceremonial heart of the PCF‑SEA forge.

---
