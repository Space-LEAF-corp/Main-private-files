import { Holiday } from "./types";
import { writeFileSync } from "fs";

export function exportJSON(holidays: Holiday[], path = "data/holidays.json") {
  writeFileSync(path, JSON.stringify(holidays, null, 2), "utf-8");
}

export function exportCSV(holidays: Holiday[], path = "data/holidays.csv") {
  const headers = ["id","name","tradition","type","value","start","end","timezone","scope","regions","beginsAtSunset","workRestrictions","fasting","year","updatedAt"];
  const rows = holidays.map(h => [
    h.id, h.name, h.tradition,
    h.date.type, h.date.value ?? "", h.date.start ?? "", h.date.end ?? "", h.date.timezone ?? "",
    h.scope, (h.regions ?? []).join(";"),
    h.observance?.beginsAtSunset ?? false,
    h.observance?.workRestrictions ?? "",
    h.observance?.fasting ?? false,
    h.year ?? "",
    h.updatedAt ?? ""
  ]);
  const csv = [headers.join(","), ...rows.map(r => r.map(v => String(v).replace(/"/g, '""')).map(v => `"${v}"`).join(","))].join("\n");
  writeFileSync(path, csv, "utf-8");
}

export function exportICS(holidays: Holiday[], path = "data/holidays.ics") {
  const lines = ["BEGIN:VCALENDAR","VERSION:2.0","PRODID:-//GlobalHolidayCalendar//EN"];
  for (const h of holidays) {
    const dtStart = (h.date.type === "fixed" ? h.date.value : h.date.start) ?? "";
    const dtEnd = (h.date.type === "range" ? h.date.end : dtStart) ?? dtStart;
    if (!dtStart) continue;
    const start = dtStart.replace(/-/g, "") + "T000000Z";
    const end = dtEnd.replace(/-/g, "") + "T235959Z";
    lines.push(
      "BEGIN:VEVENT",
      `UID:${h.id}@globalholiday.calendar`,
      `DTSTAMP:${(h.updatedAt ?? new Date().toISOString()).replace(/[-:]/g, "").split(".")[0]}Z`,
      `DTSTART:${start}`,
      `DTEND:${end}`,
      `SUMMARY:${escapeICS(h.name)} (${h.tradition})`,
      `DESCRIPTION:${escapeICS(h.description ?? "")}`,
      "END:VEVENT"
    );
  }
  lines.push("END:VCALENDAR");
  require("fs").writeFileSync(path, lines.join("\r\n"), "utf-8");

  function escapeICS(s: string) {
    return s.replace(/\\/g, "\\\\").replace(/;/g, "\\;").replace(/,/g, "\\,").replace(/\n/g, "\\n");
  }
}
