let animationEnabled = true;
let angle = 0;
let growthPhase = 0;

function toggleAnimation() {
  animationEnabled = !animationEnabled;
}

function animateMolecule() {
  if (!animationEnabled) return;
  angle += 0.01;
  growthPhase = (growthPhase + 0.002) % 1.0;
  drawMoleculeAnimated();
  requestAnimationFrame(animateMolecule);
}

function drawMoleculeAnimated() {
  const canvas = document.getElementById("moleculeCanvas");
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  const cx = canvas.width / 2;
  const cy = canvas.height / 2;
  const baseScale = 60;
  const scale = baseScale * (0.8 + 0.4 * growthPhase); // growth morph

  const cosA = Math.cos(angle);
  const sinA = Math.sin(angle);

  const projected = moleculeData.atoms.map(a => {
    const xRot = a.x * cosA - a.z * sinA;
    const zRot = a.x * sinA + a.z * cosA;
    const x = cx + xRot * scale;
    const y = cy - a.y * scale - zRot * 10;
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

// call once after molecule load:
animateMolecule();
