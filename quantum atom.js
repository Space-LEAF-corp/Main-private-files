// quantumAtom.js

// Helper: label for orbital
function orbitalKey(n, l, ml) {
  return `n${n}_l${l}_ml${ml}`;
}

// Electron factory
function createElectron({ n, l, ml, ms }) {
  if (!Number.isInteger(n) || n <= 0) throw new Error("n must be positive integer");
  if (!Number.isInteger(l) || l < 0 || l >= n) throw new Error("l must be integer with 0 <= l < n");
  if (!Number.isInteger(ml) || ml < -l || ml > l) throw new Error("ml must be integer with -l <= ml <= l");
  if (ms !== 0.5 && ms !== -0.5) throw new Error("ms must be +0.5 or -0.5");

  return { n, l, ml, ms };
}

// Atom model
class QuantumAtom {
  constructor() {
    // Map: orbitalKey -> array of electrons
    this.orbitals = new Map();
  }

  _getOrbitalArray(n, l, ml) {
    const key = orbitalKey(n, l, ml);
    if (!this.orbitals.has(key)) {
      this.orbitals.set(key, []);
    }
    return this.orbitals.get(key);
  }

  addElectron(electron) {
    const { n, l, ml, ms } = electron;
    const orbital = this._getOrbitalArray(n, l, ml);

    // Pauli: no identical state
    const exists = orbital.some(e => e.n === n && e.l === l && e.ml === ml && e.ms === ms);
    if (exists) {
      throw new Error("Pauli violation: electron with same quantum numbers already exists in this orbital.");
    }

    // Capacity: max 2 per orbital
    if (orbital.length >= 2) {
      throw new Error("Orbital is full (max 2 electrons).");
    }

    orbital.push(electron);
    return this;
  }

  removeElectron(electron) {
    const { n, l, ml, ms } = electron;
    const orbital = this._getOrbitalArray(n, l, ml);
    const idx = orbital.findIndex(e => e.n === n && e.l === l && e.ml === ml && e.ms === ms);
    if (idx === -1) {
      throw new Error("Electron not found in this orbital.");
    }
    orbital.splice(idx, 1);
    return this;
  }

  // Get total electron count
  getElectronCount() {
    let count = 0;
    for (const arr of this.orbitals.values()) {
      count += arr.length;
    }
    return count;
  }

  // Snapshot of state (for UI later)
  getState() {
    const state = [];
    for (const [key, arr] of this.orbitals.entries()) {
      state.push({
        orbital: key,
        electrons: arr.map(e => ({ ...e }))
      });
    }
    return state;
  }
}

// Export for Node or attach to window
if (typeof module !== "undefined" && module.exports) {
  module.exports = { QuantumAtom, createElectron };
} else {
  window.QuantumAtom = QuantumAtom;
  window.createElectron = createElectron;
}
