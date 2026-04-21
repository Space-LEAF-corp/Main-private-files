// tests.js
// If in Node:
// const { QuantumAtom, createElectron } = require("./quantumAtom");

function runTests() {
  const atom = new QuantumAtom();

  // Test 1: add two electrons to same orbital with opposite spin
  const e1 = createElectron({ n: 1, l: 0, ml: 0, ms: 0.5 });
  const e2 = createElectron({ n: 1, l: 0, ml: 0, ms: -0.5 });

  atom.addElectron(e1);
  atom.addElectron(e2);

  console.log("Test 1 passed:", atom.getElectronCount() === 2);

  // Test 2: Pauli violation (same ms)
  let pauliErrorCaught = false;
  try {
    const e3 = createElectron({ n: 1, l: 0, ml: 0, ms: 0.5 });
    atom.addElectron(e3);
  } catch (err) {
    pauliErrorCaught = true;
    console.log("Test 2 (Pauli) error message:", err.message);
  }
  console.log("Test 2 passed:", pauliErrorCaught);

  // Test 3: orbital capacity (third electron in same orbital)
  let capacityErrorCaught = false;
  try {
    const atom2 = new QuantumAtom();
    atom2.addElectron(e1);
    atom2.addElectron(e2);
    const e4 = createElectron({ n: 1, l: 0, ml: 0, ms: -0.5 });
    atom2.addElectron(e4);
  } catch (err) {
    capacityErrorCaught = true;
    console.log("Test 3 (capacity) error message:", err.message);
  }
  console.log("Test 3 passed:", capacityErrorCaught);

  // Test 4: different orbitals
  const atom3 = new QuantumAtom();
  atom3.addElectron(createElectron({ n: 2, l: 1, ml: -1, ms: 0.5 }));
  atom3.addElectron(createElectron({ n: 2, l: 1, ml: 0, ms: 0.5 }));
  atom3.addElectron(createElectron({ n: 2, l: 1, ml: 1, ms: -0.5 }));

  console.log("Test 4 passed:", atom3.getElectronCount() === 3);

  console.log("Final state example:", JSON.stringify(atom3.getState(), null, 2));
}

runTests();
