const seedLibrary = {
  "Class I — Ocean-Born Hybrid": {
    molecule_name: "Hybrid Seed I",
    atoms: [
      { element: "C", x: 0, y: 0, z: 0 },
      { element: "O", x: 1, y: 0, z: 0 },
      { element: "Si", x: -1, y: 0, z: 0 }
    ],
    bonds: [
      { a: 0, b: 1, order: 1 },
      { a: 0, b: 2, order: 1 }
    ]
  },
  "Class II — High-Pressure Inorganic": {
    molecule_name: "Inorganic Seed II",
    atoms: [
      { element: "Si", x: 0, y: 0, z: 0 },
      { element: "O", x: 1, y: 1, z: 0 },
      { element: "O", x: -1, y: -1, z: 0 },
      { element: "O", x: 1, y: -1, z: 0 },
      { element: "O", x: -1, y: 1, z: 0 }
    ],
    bonds: [
      { a: 0, b: 1, order: 1 },
      { a: 0, b: 2, order: 1 },
      { a: 0, b: 3, order: 1 },
      { a: 0, b: 4, order: 1 }
    ]
  },
  "Class III — Lineage / Premium Core": {
    molecule_name: "Lineage Seed III",
    atoms: [
      { element: "C", x: 0, y: 0, z: 0 },
      { element: "C", x: 1, y: 1, z: 1 },
      { element: "C", x: -1, y: -1, z: -1 },
      { element: "O", x: 2, y: 0, z: 0 }
    ],
    bonds: [
      { a: 0, b: 1, order: 2 },
      { a: 0, b: 2, order: 2 },
      { a: 0, b: 3, order: 1 }
    ]
  }
};

function loadSeedForClass(className) {
  const seed = seedLibrary[className];
  if (!seed) return;
  moleculeData = { atoms: seed.atoms, bonds: seed.bonds };
  renderAtoms();
  renderBonds();
  drawMoleculeAnimated();
}
