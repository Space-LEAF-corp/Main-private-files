# Ceremonial Terrain v1 – Map Generator Prototype

Ceremonial Terrain v1 is a prototype **random map generator** and **planet simulator** designed for
kids, game designers, and stewards. It turns abstract patterns (like tea or coffee residue, noise,
or hand‑drawn maps) into:

- 3D terrain (surface + ocean floor)
- 4D evolution over time (tectonic plates + volcanic activity)
- Educational analysis views (plates, volcanoes, bathymetry)

It also includes **application‑level security scaffolding** to help protect users' projects.

---

## Features

- **Terrain engine**
  - Heightmap‑based surface generation
  - Ocean floor / bathymetry layer
  - Tectonic plate partitioning and motion
  - Volcanic hotspot and boundary‑driven eruptions

- **Simulation**
  - Deterministic 365‑step (day/epoch) world evolution
  - Forward and backward time stepping using snapshots
  - Basic validation and consistency checks

- **Security scaffolding**
  - Simple in‑memory user accounts
  - Project integrity hashing (tamper detection)
  - Hooks for encryption and stricter auth

- **Extensibility**
  - Hooks for image‑based input (tea leaves, drawings)
  - Hooks for 3D/4D renderers (e.g., game engines, WebGL, Blender)

> NOTE: This is a **prototype** meant as a reference/starting point. It is **not** a complete
security solution or production‑grade engine.

---

## Installation

```bash
git clone https://github.com/<your-username>/ceremonial-terrain-v1.git
cd ceremonial-terrain-v1
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -e .
