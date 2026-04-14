// server.js
import express from "express";
import crypto from "crypto";

const app = express();
app.use(express.json());

const elections = new Map();
const ballots = new Map();
const audits = new Map();

function sha256(obj) {
  return crypto.createHash("sha256")
    .update(JSON.stringify(obj))
    .digest("hex");
}

function appendAudit(electionId, type, payload) {
  const chain = audits.get(electionId) || [];
  const prevHash = chain.length ? chain[chain.length - 1].entryHash : null;
  const entry = {
    id: crypto.randomUUID(),
    electionId,
    type,
    timestamp: new Date().toISOString(),
    payloadHash: sha256(payload),
    prevEntryHash: prevHash
  };
  entry.entryHash = sha256(entry);
  chain.push(entry);
  audits.set(electionId, chain);
}

app.post("/api/elections", (req, res) => {
  const id = req.body.id;
  if (!id) return res.status(400).json({ error: "id required" });
  elections.set(id, {
    id,
    name: req.body.name,
    status: "open",
    createdAt: new Date().toISOString(),
    seal: null
  });
  ballots.set(id, []);
  appendAudit(id, "ELECTION_CREATED", { id, name: req.body.name });
  res.status(201).json({ id });
});

app.post("/api/elections/:id/ballots", (req, res) => {
  const election = elections.get(req.params.id);
  if (!election || election.status !== "open") {
    return res.status(400).json({ error: "Election not open" });
  }
  const ballot = {
    id: crypto.randomUUID(),
    electionId: election.id,
    choices: req.body.choices,
    castAt: new Date().toISOString(),
    districtId: req.body.districtId,
    method: req.body.method,
    demographicBucket: req.body.demographicBucket
  };
  ballot.hash = sha256(ballot);
  ballots.get(election.id).push(ballot);
  appendAudit(election.id, "BALLOT_CAST", { ballotHash: ballot.hash });
  res.status(201).json({ ballotId: ballot.id });
});

app.post("/api/elections/:id/close", (req, res) => {
  const election = elections.get(req.params.id);
  if (!election) return res.status(404).json({ error: "Not found" });

  election.status = "closed";

  const b = ballots.get(election.id) || [];

  // Aggregate turnout streams
  const districts = {};
  const methodsAgg = {};
  const demoAgg = {};

  for (const ballot of b) {
    const d = ballot.districtId;
    const m = ballot.method;
    const g = ballot.demographicBucket;

    districts[d] = districts[d] || { districtId: d, registered: 0, ballotsCast: 0 };
    districts[d].ballotsCast++;

    methodsAgg[m] = methodsAgg[m] || { method: m, registered: 0, ballotsCast: 0 };
    methodsAgg[m].ballotsCast++;

    demoAgg[g] = demoAgg[g] || { group: g, registered: 0, ballotsCast: 0 };
    demoAgg[g].ballotsCast++;
  }

  const results = {
    districts: Object.values(districts),
    methods: Object.values(methodsAgg),
    demographics: Object.values(demoAgg)
  };

  const seal = {
    schemaVersion: "1.0",
    resultsHash: sha256(results),
    ballotSetHash: sha256(b.map(x => x.hash)),
    signature: null // plug in real signing later
  };

  election.seal = seal;
  appendAudit(election.id, "SEAL_ISSUED", seal);

  res.json({ election, results });
});

app.get("/api/elections/:id/results", (req, res) => {
  const election = elections.get(req.params.id);
  if (!election || !election.seal) {
    return res.status(404).json({ error: "No results yet" });
  }
  const b = ballots.get(election.id) || [];
  // recompute aggregates as above or store them
  // (for brevity, assume we stored them in election.results)
  res.json({
    election: {
      id: election.id,
      name: election.name,
      seal: election.seal
    },
    results: election.results
  });
});

app.get("/api/elections/:id/audit", (req, res) => {
  res.json(audits.get(req.params.id) || []);
});

app.listen(3000, () => console.log("Clear Voting backend on :3000"));
