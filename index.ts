import { validateHoliday } from "./core/validate";
import { normalize } from "./core/normalize";
import { mergeHolidays } from "./core/merge";
import { exportJSON, exportCSV, exportICS } from "./core/exports";
import { adapters } from "./adapters/user-contrib"; // plus static adapters
import { staticJewish } from "./adapters/static-jewish";
import { staticIslamic } from "./adapters/static-islamic";
import { staticChristian } from "./adapters/static-christian";
import { staticBahai } from "./adapters/static-bahai";
import { staticHindu } from "./adapters/static-hindu";
import { staticBuddhist } from "./adapters/static-buddhist";

async function build() {
  const allSources = [];
  // Run adapters concurrently
  const results = await Promise.all([
    staticJewish(), staticIslamic(), staticChristian(),
    staticBahai(), staticHindu(), staticBuddhist(),
    ...adapters.map(a => a.run())
  ]);

  for (const res of results) {
    const cleaned = res.map(h => normalize(validateHoliday(h)));
    allSources.push(cleaned);
  }

  const merged = mergeHolidays(allSources);
  exportJSON(merged);
  exportCSV(merged);
  exportICS(merged);
  console.log(`Built ${merged.length} holidays -> data/holidays.*`);
}

build().catch(err => {
  console.error(err);
  process.exit(1);
});
