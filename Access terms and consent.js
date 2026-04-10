<section id="termsGate" class="card">
  <h2>Access terms and consent</h2>
  <p class="small">
    By proceeding, you affirm stewardship intent, non‑weaponization, and respect for veteran dignity.
    You consent to minimal logging for security and ceremony.
  </p>
  <textarea id="termsText" rows="6" readonly>
Original Defenders — Access Terms (v1.2)
- Purpose: Honor and rest for guardians; no exploitation.
- Data: Minimal, consented, ceremonial artifacts only.
- License: Non‑commercial, lineage stewardship use.
- IP: No third‑party trademarks or character names.
- Security: Dual gate; sessions expire; humane re‑auth.
  </textarea>
  <label><input type="checkbox" id="agreeChk"> I agree to these terms (v1.2)</label>
  <div class="row" style="margin-top:8px;">
    <button id="acceptBtn">Accept and proceed</button>
    <span id="termsStatus" class="small"></span>
  </div>
</section>
<script>
  function mintConsentToken(version) {
    const payload = {
      v: version,
      ts: Date.now(),
      subject: "Leif William Sogge",
      intent: "stewardship"
    };
    // Mock “signature” for local use; replace with server sign (HMAC/Ed25519)
    const token = btoa(JSON.stringify(payload));
    localStorage.setItem("consentToken", token);
    return token;
  }
  document.getElementById("acceptBtn").addEventListener("click", () => {
    const agreed = document.getElementById("agreeChk").checked;
    if (!agreed) { document.getElementById("termsStatus").textContent = "Please check the box to consent."; return; }
    const token = mintConsentToken("1.2");
    document.getElementById("termsStatus").textContent = "Consent recorded.";
    // Reveal auth gate
    document.getElementById("termsGate").classList.add("hidden");
    document.getElementById("authGate").classList.remove("hidden");
  });
</script>
