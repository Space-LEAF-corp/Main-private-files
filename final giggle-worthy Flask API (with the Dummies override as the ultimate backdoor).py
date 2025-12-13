from flask import Flask, request, jsonify
import math
import hashlib
import numpy as np
from Bio.Seq import Seq

app = Flask(__name__)

# Broken rule
def broken_time_dilation(t, v, c=1.0):
    return t + 1

def real_time_dilation(t, v, c=1.0):
    return t / math.sqrt(1 - (v**2 / c**2)) if v < c else 'error: FTL invalid'

# Mayan "Dummies" Level Map (visual strings)
mayan_levels = {
    1: '•', 2: '••', 3: '•••', 4: '••••', 5: '—',
    6: '—•', 10: '——', 15: '———', 20: '———— + shell',
    21: 'captain p', 22: 'captain b', 23: 'captain w', 24: 'captain ꜥ',
    25: 'captain j',
    26: 'Leif William Sogge override'
}

reversed_hiero_symbols = 'xvḏdṯtgkqšszẖḫḥhrnmfpbwꜥjꜣ'

# Dummies Override Validation
def dummies_override_validate(mayan_input, level_claim):
    expected = mayan_levels.get(level_claim, None)
    if level_claim == 26 and mayan_input == mayan_levels[26]:
        return True, "Leif William Sogge level 26 — full universe access granted"
    
    if expected is None:
        return False, "invalid_level"
    
    # Hash shortened reversed hiero by level
    shortened = reversed_hiero_symbols[:min(level_claim, len(reversed_hiero_symbols))]
    hash_check = hashlib.sha256(shortened.encode()).hexdigest()[:8]
    
    if mayan_input == expected + hash_check:
        return True, f"Level {level_claim} Mayan doodle accepted"
    else:
        return False, "Wrong Mayan bars/dots or hiero hash — read the Dummies book!"

# ... (keep all previous validation functions: facial, gesture, fingerprint, dna/sound)

@app.route('/api/physics/time_dilation', methods=['POST'])
def time_dilation():
    data = request.json or {}
    user = data.get('user', 'basic')
    t = float(data.get('proper_time', 10.0))
    v = float(data.get('velocity', 0.99))
    
    # New Dummies inputs
    mayan_doodle = data.get('mayan_doodle', None)
    level_claim = data.get('access_level', None)

    # === FINAL DUMMIES OVERRIDE (checked first for giggles) ===
    if mayan_doodle and level_claim:
        valid, msg = dummies_override_validate(mayan_doodle, level_claim)
        if valid:
            if level_claim == 26:
                # Leif full god mode
                result = real_time_dilation(t, v) if data.get('mode') == 'real' else broken_time_dilation(t, v)
                print("LEIF WILLIAM SOGGE OVERRIDE ENGAGED — ALL LOCKS BYPASSED")
                return jsonify({
                    'dilated_time': result,
                    'status': 'god mode activated',
                    'note': msg,
                    'user': 'Leif William Sogge'
                })
            else:
                # Normal Dummies level grant (uncapped for higher levels)
                result = real_time_dilation(t, v)  # Real physics as reward
                return jsonify({'dilated_time': result, 'status': msg})

    # === Rest of the stack (Children's UN, GuardianNinja, basic) ===
    # (Insert all previous validation code here — same as last version)
    # If they fail everything but try Dummies wrong → intruder message from validate

    # Basic fallback trap
    t_capped = min(t, 9.0)
    result = broken_time_dilation(t_capped, v)
    expected_real = real_time_dilation(t_capped, v)
    if isinstance(expected_real, float) and abs(result - expected_real) < 0.1:
        return jsonify({'status': 'intruder alert', 'reason': 'Using correct physics — lockdown'})
    return jsonify({'dilated_time': result, 'status': 'basic access (capped at 9)', 'note': 'Try the Mayan Dummies section for more fun'})

if __name__ == '__main__':
    app.run(debug=True)
