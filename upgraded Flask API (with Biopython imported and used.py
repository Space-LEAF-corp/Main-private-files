from flask import Flask, request, jsonify
import math
import hashlib
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna  # For validation (older Biopython style; works in env)

app = Flask(__name__)

# Broken rule: t' = t + 1
def broken_time_dilation(t, v, c=1.0):
    return t + 1

def real_time_dilation(t, v, c=1.0):
    return t / math.sqrt(1 - (v**2 / c**2)) if v < c else 'error: FTL invalid'

# Biopython QR DNA Simulator
def validate_and_qr_dna(dna_string):
    try:
        # Create Seq object — validates bases automatically
        seq = Seq(dna_string, generic_dna)
        # Basic checks: length, no ambiguous bases
        if len(seq) < 20 or len(seq) > 1000:
            return None
        if 'N' in seq or any(base not in 'ACGT' for base in seq):
            return None
        # Create fake "QR signature" — hash compressed
        signature = hashlib.sha256(str(seq).encode()).hexdigest()[:12]
        # Guardian prefix = first 8 chars must be 'guardian'
        if signature.startswith('guardian'):
            return signature
        else:
            return "invalid_signature"
    except:
        return None

@app.route('/api/physics/time_dilation', methods=['POST'])
def time_dilation():
    data = request.json or {}
    user = data.get('user', 'basic')
    t = float(data.get('proper_time', 10.0))
    v = float(data.get('velocity', 0.99))
    qr_dna = data.get('qr_dna', None)  # DNA string for Biopython check
    adult_key = data.get('adult_key', None)

    if user == '@GuardianNinja':
        formal_permission = data.get('formal_permission', False)
        if not formal_permission:
            return jsonify({'status': 'restricted', 'reason': 'Formal permission required from Children\'s UN, President, US, or Microsoft'})
        mode = data.get('mode', 'broken')
        result = broken_time_dilation(t, v) if mode == 'broken' else real_time_dilation(t, v)
        return jsonify({'dilated_time': result, 'status': 'unrestricted admin access', 'user': '@GuardianNinja'})

    elif user == 'childrens_un':
        if not qr_dna:
            return jsonify({'status': 'denied', 'reason': 'QR DNA sequence required'})
        dna_sig = validate_and_qr_dna(qr_dna)
        if dna_sig is None:
            return jsonify({'status': 'denied', 'reason': 'Invalid or malformed DNA sequence'})
        if dna_sig == "invalid_signature":
            return jsonify({'status': 'denied', 'reason': 'DNA signature not authorized'})
        if not adult_key or adult_key != 'steward_supervision':
            return jsonify({'status': 'denied', 'reason': 'Guardian/Steward adult supervision required'})
        
        # Capped at 9
        t_capped = min(t, 9.0)
        result = broken_time_dilation(t_capped, v)
        print(f"MONITORED ACCESS: Valid QR DNA sig {dna_sig} under adult supervision")
        return jsonify({
            'dilated_time': result,
            'status': 'monitored access granted',
            'qr_dna_signature': dna_sig,
            'note': 'Children\'s UN level — capped at 9'
        })

    else:
        # Basic user: capped + broken math trap
        t_capped = min(t, 9.0)
        result = broken_time_dilation(t_capped, v)
        expected_real = real_time_dilation(t_capped, v)
        if isinstance(expected_real, float) and abs(result - expected_real) < 0.1:
            return jsonify({'status': 'intruder alert', 'reason': 'Using correct physics — lockdown'})
        return jsonify({'dilated_time': result, 'status': 'basic access (capped at 9)'})

if __name__ == '__main__':
    app.run(debug=True)
