import { Adapter } from "./adapter-template";
import { readFileSync, existsSync } from "fs";
import path from "path";

export const adapters: Adapter[] = [
  {
    name: "user-contrib-json",
    run: async () => {
      const p = path.resolve("contrib/examples/sample-contrib.json");
      if (!existsSync(p)) return [];
      const raw = JSON.parse(readFileSync(p, "utf-8"));
      return Array.isArray(raw) ? raw : [];
    }
  }
];
