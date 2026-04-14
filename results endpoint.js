async function loadResults(electionId) {
  const res = await fetch(`/api/elections/${electionId}/results`);
  if (!res.ok) throw new Error("Failed to load results");
  return res.json();
}

async function init() {
  const { election, results } = await loadResults("2026-ocala-municipal");
  const { districts, methods, demographics } = results;

  // existing sort controls...
  // then renderTable(...) as you already do

  // show seal in header/footer
  const lineageEl = document.querySelector(".lineage");
  lineageEl.textContent =
    `Authorship: Leif William Sogge • Version ${election.seal.schemaVersion} • ` +
    `Results hash: ${election.seal.resultsHash.slice(0, 12)}…`;
}
