import { runAllSimulations } from "./simulation";

(async () => {
  try {
    const results = await runAllSimulations();
    console.log("Simulation Results:", results);
  } catch (err) {
    console.error("Simulation Error:", err);
    process.exit(1);
  }
})();
