{
  "artifact": {
    "name": "Clear Voting Seal",
    "version": "1.0.0",
    "authors": [
      {
        "name": "Leif William Sogge",
        "role": "Ceremonial steward",
        "provenance": "Ark Captain’s Log"
      }
    ],
    "purpose": "Separate turnout rates into distinct, transparent streams with unified lineage.",
    "lineage": {
      "parent": "Ark Governance Framework",
      "anchors": [
        "Resilience under interruption",
        "Explicit authorship",
        "Community visibility"
      ]
    },
    "integrity": {
      "hashing": "SHA-256",
      "signing": "Detached signature optional",
      "schemaVersion": "1.0"
    }
  },
  "display": {
    "sections": ["Districts", "Methods", "Demographics", "Lineage"],
    "defaults": {
      "sortBy": "turnoutRate",
      "order": "desc",
      "thresholds": {
        "low": 0.35,
        "medium": 0.55,
        "high": 0.70
      }
    }
  },
  "dataRefs": {
    "districts": "./data/districts.json",
    "methods": "./data/methods.json",
    "demographics": "./data/demographics.json"
  }
}
[
  { "districtId": "FL-05-001", "name": "Ocala Precinct 1", "registered": 1200, "ballotsCast": 804 },
  { "districtId": "FL-05-002", "name": "Ocala Precinct 2", "registered": 1475, "ballotsCast": 876 },
  { "districtId": "FL-05-003", "name": "Marion Rural A", "registered": 980, "ballotsCast": 522 }
]
[
  { "group": "18-29", "registered": 900, "ballotsCast": 468 },
  { "group": "30-44", "registered": 1100, "ballotsCast": 638 },
  { "group": "45-64", "registered": 980, "ballotsCast": 606 },
  { "group": "65+", "registered": 675, "ballotsCast": 484 }
]
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Clear Voting Seal</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <style>
    :root { --low:#ffb3b3; --med:#ffd67a; --high:#a6f3a6; --ink:#1c1c1c; --bg:#f7f7f9; }
    body { font-family: system-ui, -apple-system, Segoe UI, Roboto, sans-serif; background: var(--bg); color: var(--ink); margin: 0; }
    header { padding: 16px 20px; background: #ffffff; border-bottom: 1px solid #ddd; }
    h1 { margin: 0 0 4px 0; font-size: 20px; }
    .lineage { font-size: 12px; color: #555; }
    .grid { display: grid; grid-template-columns: repeat(auto-fit,minmax(280px,1fr)); gap: 16px; padding: 16px; }
    .card { background: #fff; border: 1px solid #e6e6e6; border-radius: 10px; overflow: hidden; }
    .card h2 { font-size: 16px; margin: 0; padding: 12px 12px 0; }
    .table { width: 100%; border-collapse: collapse; }
    .table th, .table td { text-align: left; padding: 10px 12px; border-top: 1px solid #eee; font-size: 14px; }
    .rate { display: inline-block; padding: 4px 8px; border-radius: 6px; font-weight: 600; }
    .low { background: var(--low); }
    .med { background: var(--med); }
    .high { background: var(--high); }
    footer { padding: 16px 20px; font-size: 12px; color: #666; }
    .badge { display:inline-block; margin-right:8px; padding:4px 6px; border:1px solid #ddd; border-radius:6px; font-size:12px; }
    .controls { padding: 0 16px 16px; display:flex; gap:8px; align-items:center; }
    select { padding:6px 8px; }
  </style>
</head>
<body>
  <header>
    <h1>Clear Voting Seal — separated rates with unified lineage</h1>
    <div class="lineage">Authorship: Leif William Sogge • Version 1.0.0 • Ark Governance Framework</div>
    <div class="controls">
      <label>Sort by:
        <select id="sort">
          <option value="rate">Turnout rate</option>
          <option value="cast">Ballots cast</option>
          <option value="registered">Registered</option>
        </select>
      </label>
      <label>Order:
        <select id="order">
          <option value="desc">Descending</option>
          <option value="asc">Ascending</option>
        </select>
      </label>
    </div>
    <div>
      <span class="badge">Low &lt; 35%</span>
      <span class="badge">Medium 35–55%</span>
      <span class="badge">High &gt; 55%</span>
    </div>
  </header>

  <main class="grid">
    <section class="card">
      <h2>Districts</h2>
      <table class="table" id="districts"></table>
    </section>

    <section class="card">
      <h2>Methods</h2>
      <table class="table" id="methods"></table>
    </section>

    <section class="card">
      <h2>Demographics</h2>
      <table class="table" id="demographics"></table>
    </section>
  </main>

  <footer>
    Unified context: these separated rates form one election lineage. Separation is for clarity, not division.
  </footer>

  <script>
    const thresholds = { low: 0.35, med: 0.55 };
    async function loadJSON(path) {
      const res = await fetch(path);
      if (!res.ok) throw new Error("Failed to load " + path);
      return res.json();
    }
    function rateClass(r) {
      if (r < thresholds.low) return "low";
      if (r < thresholds.med) return "med";
      return "high";
    }
    function fmtRate(r) { return (r * 100).toFixed(1) + "%"; }
    function sortData(arr, by, order) {
      const key = by === "rate" ? (o => (o.ballotsCast / o.registered)) :
                  by === "cast" ? (o => o.ballotsCast) :
                                  (o => o.registered);
      const s = arr.slice().sort((a,b) => key(b) - key(a));
      return order === "asc" ? s.reverse() : s;
    }
    function renderTable(el, rows, columns, sortBy, order) {
      const sorted = sortData(rows, sortBy, order);
      el.innerHTML = `
        <tr>
          ${columns.map(c => `<th>${c.label}</th>`).join("")}
          <th>Turnout</th>
        </tr>
        ${sorted.map(row => {
          const rate = row.ballotsCast / row.registered;
          return `<tr>
            ${columns.map(c => `<td>${row[c.key]}</td>`).join("")}
            <td><span class="rate ${rateClass(rate)}">${fmtRate(rate)}</span></td>
          </tr>`;
        }).join("")}
      `;
    }
    async function init() {
      const [districts, methods, demographics] = await Promise.all([
        loadJSON("./data/districts.json"),
        loadJSON("./data/methods.json"),
        loadJSON("./data/demographics.json")
      ]);
      const sortSel = document.getElementById("sort");
      const orderSel = document.getElementById("order");
      function draw() {
        const sortBy = sortSel.value;
        const order = orderSel.value;
        renderTable(document.getElementById("districts"), districts,
                    [{label:"District", key:"name"}, {label:"Registered", key:"registered"}, {label:"Cast", key:"ballotsCast"}],
                    sortBy, order);
        renderTable(document.getElementById("methods"), methods,
                    [{label:"Method", key:"method"}, {label:"Registered", key:"registered"}, {label:"Cast", key:"ballotsCast"}],
                    sortBy, order);
        renderTable(document.getElementById("demographics"), demographics,
                    [{label:"Group", key:"group"}, {label:"Registered", key:"registered"}, {label:"Cast", key:"ballotsCast"}],
                    sortBy, order);
      }
      sortSel.addEventListener("change", draw);
      orderSel.addEventListener("change", draw);
      draw();
    }
    init().catch(err => {
      console.error(err);
      alert("Error loading data. Check file paths.");
    });
  </script>
</body>
</html>
