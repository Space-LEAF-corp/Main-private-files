import { Holiday } from "../core/types";

export interface Adapter {
  name: string;
  run: () => Promise<Holiday[]>;
}

export function createStaticAdapter(name: string, holidays: Holiday[]): Adapter {
  return {
    name,
    run: async () => holidays
  };
}


Static examples (e.g., src/adapters/static-jewish.ts)

import { createStaticAdapter } from "./adapter-template";
import { Holiday } from "../core/types";

const year = new Date().getUTCFullYear();

const holidays: Holiday[] = [
  {
    id: `shabbat-${year}-weekly`,
    name: "Shabbat",
    aliases: ["Sabbath"],
    tradition: "Jewish",
    date: { type: "rule", value: "WEEKLY-FRIDAY-SUNSET-TO-SATURDAY-NIGHTFALL" },
    observance: { beginsAtSunset: true, workRestrictions: "Traditional restrictions apply" },
    scope: "global",
    description: "Weekly day of rest and holiness from Friday evening to Saturday night.",
    year
  }
  // Add computed date adapters later for specific annual holidays (e.g., Rosh Hashanah)
];

export async function staticJewish() {
  return createStaticAdapter("static-jewish", holidays).run();
}
