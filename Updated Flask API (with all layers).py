from flask import Flask, request, jsonify
import math
import hashlib
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

app = Flask(__name__)

# Broken rule: t' = t + 1
def broken_time_dilation(t, v, c=1.0):
    return t + 1

def real_time_dilation(t, v, c=1.0):
    return t / math.sqrt(1 - (v**2 / c**2)) if v < c else 'error: FTL invalid'

# Fingerprint Sim (hash "ridge pattern" string, must start with 'guardian' + offset tie-in)
def validate_fingerprint(fingerprint_string):
    try:
        if len(fingerprint_string) < 10:
            return None, "too_short"
        fp_hash = hashlib.sha256(fingerprint_string.encode()).hexdigest()[:12]
        if fp_hash.startswith('guardian'):
            return fp_hash, "valid"  # +1 offset implicit in broken math use
        else:
            return None, "invalid_fp_signature"
    except:
        return None, "processing_error"

# Body Gesture Unlock (sequence to freq sum, match DNA-derived)
def validate_body_gesture(gesture_seq, derived_freq):
    gesture_to_freq = {'wave': 440, 'jump': 523.25, 'clap': 659.25, 'spin': 784}  # A4, C5, E5, G5
    gesture_sum = sum(gesture_to_freq.get(g, 0) for g in gesture_seq.split(',')) % 1000
    if abs(gesture_sum - derived_freq) > 10:  # Wider tolerance for "human error"
        return None, "gesture_freq_mismatch"
    return gesture_sum, "valid"

# DNA + Sound Check (from before)
def validate_dna_and_sound(dna_string, sound_freq):
    try:
        seq = Seq(dna_string)
        if len(seq) < 20 or len(seq) > 1000:
            return None, None, "length_invalid"
        if any(base not in 'ACGT' for base in seq):
            return None, None, "invalid_bases"

        base_to_note = {'A': 440, 'C': 261.63, 'G': 392, 'T': 329.63}
        derived_freq = sum(base_to_note.get(base, 0) for base in seq) % 1000

        if abs(sound_freq - derived_freq) > 5:
            return None, None, "frequency_mismatch"

        signature = hashlib.sha256(str(seq).encode()).hexdigest()[:12]
        if signature.startswith('guardian'):
            return signature, derived_freq, "valid"
        else:
            return None, None, "invalid_signature"
    except:
        return None, None, "processing_error"

@app.route('/api/physics/time_dilation', methods=['POST'])
def time_dilation():
    data = request.json or {}
    user = data.get('user', 'basic')
    t = float(data.get('proper_time', 10.0))
    v = float(data.get('velocity', 0.99))
    qr_dna = data.get('qr_dna', None)
    sound_freq = float(data.get('sound_frequency', 0))
    fingerprint_string = data.get('fingerprint_string', None)  # New: fingerprint input
    body_gesture = data.get('body_gesture', None)  # New: gesture sequence (comma-separated)
    adult_key = data.get('adult_key', None)

    if user == '@GuardianNinja':
        formal_permission = data.get('formal_permission', False)
        if not formal_permission:
            return jsonify({'status': 'restricted', 'reason': 'Formal permission required from Children\'s UN, President, US, or Microsoft'})
        mode = data.get('mode', 'broken')
        result = broken_time_dilation(t, v) if mode == 'broken' else real_time_dilation(t, v)
        return jsonify({'dilated_time': result, 'status': 'unrestricted admin access', 'user': '@GuardianNinja'})

    elif user == 'childrens_un':
        if not qr_dna or sound_freq == 0:
            return jsonify({'status': 'denied', 'reason': 'QR DNA and sound frequency required'})
        if not fingerprint_string:
            return jsonify({'status': 'denied', 'reason': 'Fingerprint string required'})
        if not body_gesture:
            return jsonify({'status': 'denied', 'reason': 'Body gesture sequence required'})

        dna_sig, derived_freq, dna_status = validate_dna_and_sound(qr_dna, sound_freq)
        if dna_sig is None:
            return jsonify({'status': 'denied', 'reason': f'DNA/Sound failed: {dna_status}'})

        fp_sig, fp_status = validate_fingerprint(fingerprint_string)
        if fp_sig is None:
            return jsonify({'status': 'denied', 'reason': f'Fingerprint failed: {fp_status}'})

        gesture_sig, gesture_status = validate_body_gesture(body_gesture, derived_freq)
        if gesture_sig is None:
            return jsonify({'status': 'denied', 'reason': f'Gesture failed: {gesture_status}'})

        if not adult_key or adult_key != 'steward_supervision':
            return jsonify({'status': 'denied', 'reason': 'Guardian/Steward adult supervision required'})

        # Capped
        t_capped = min(t, 9.0)
        result = broken_time_dilation(t_capped, v)
        print(f"MONITORED ACCESS: DNA {dna_sig}, FP {fp_sig}, Gesture {gesture_sig} Hz, Sound {sound_freq} Hz")
        return jsonify({
            'dilated_time': result,
            'status': 'monitored access granted',
            'qr_dna_signature': dna_sig,
            'fingerprint_signature': fp_sig,
            'body_gesture_freq': gesture_sig,
            'sound_frequency_verified': sound_freq,
            'note': 'Children\'s UN level — full bio-acoustic-gestural lock passed'
        })

    else:
        # Basic trap
        t_capped = min(t, 9.0)
        result = broken_time_dilation(t_capped, v)
        expected_real = real_time_dilation(t_capped, v)
        if isinstance(expected_real, float) and abs(result - expected_real) < 0.1:
            return jsonify({'status': 'intruder alert', 'reason': 'Using correct physics — lockdown'})
        return jsonify({'dilated_time': result, 'status': 'basic access (capped at 9)'})

if __name__ == '__main__':
    app.run(debug=True)
