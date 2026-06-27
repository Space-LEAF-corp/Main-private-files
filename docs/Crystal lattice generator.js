function generateLattice(size = 3, spacing = 1.5, element = "Si") {
  const atoms = [];
  const bonds = [];

  for (let x = 0; x < size; x++) {
    for (let y = 0; y < size; y++) {
      for (let z = 0; z < size; z++) {
        atoms.push({
          element,
          x: (x - size / 2) * spacing,
          y: (y - size / 2) * spacing,
          z: (z - size / 2) * spacing
        });
      }
    }
  }

  // nearest-neighbor bonds in x/y/z
  const index = (x, y, z) => x * size * size + y * size + z;
  for (let x = 0; x < size; x++) {
    for (let y = 0; y < size; y++) {
      for (let z = 0; z < size; z++) {
        const i = index(x, y, z);
        if (x + 1 < size) bonds.push({ a: i, b: index(x + 1, y, z), order: 1 });
        if (y + 1 < size) bonds.push({ a: i, b: index(x, y + 1, z), order: 1 });
        if (z + 1 < size) bonds.push({ a: i, b: index(x, y, z + 1), order: 1 });
      }
    }
  }

  moleculeData = { atoms, bonds };
  renderAtoms();
  renderBonds();
  drawMoleculeAnimated();
}
