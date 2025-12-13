from flask import Flask, request, jsonify
import math
import hashlib
import numpy as np
from Bio.Seq import Seq

app = Flask(__name__)

# Broken rule: t' = t + 1
def broken_time_dilation(t, v, c=1.0):
    return t + 1

def real_time_dilation(t, v, c=1.0):
    return t / math.sqrt(1 - (v**2 / c**2)) if v < c else 'error: FTL invalid'

# Facial Recognition Sim (cosine similarity on vector)
GUARDIAN_FACE_VECTOR = np.array([0.1, 0.8, 0.3, 0.9, 0.2, 0.7, 0.5, 0.6])  # Pre-defined template

def validate_facial_recog(face_vector_list):
    try:
        face_vec = np.array(face_vector_list)
        if len(face_vec) != len(GUARDIAN_FACE_VECTOR):
            return None, "vector_length_mismatch"
        # Cosine similarity
        cosine = np.dot(face_vec, GUARDIAN_FACE_VECTOR) / (np.linalg.norm(face_vec) * np.linalg.norm(GUARDIAN_FACE_VECTOR))
        if cosine > 0.95:  # High threshold for "match"
            return cosine, "valid"
        else:
            return None, "low_similarity"
    except:
        return None, "processing_error"

# Fingerprint Sim
def validate_fingerprint(fingerprint_string):
    try:
        if len(fingerprint_string) < 10:
            return None, "too_short"
        fp_hash = hashlib.sha256(fingerprint_string.encode()).hexdigest()[:12]
        if fp_hash.startswith('guardian'):
            return fp_hash, "valid"
        else:
            return None, "invalid_fp_signature"
    except:
        return None, "processing_error"

# Body Gesture
def validate_body_gesture(gesture_seq, derived_freq):
    gesture_to_freq = {'wave': 440, 'jump': 523.25, 'clap': 659.25, 'spin': 784, 'nod': 493.88}
    gesture_sum = sum(gesture_to_freq.get(g, 0) for g in gesture_seq.split(',')) % 1000
    if abs(gesture_sum - derived_freq) > 10:
        return None, "gesture_freq_mismatch"
    return gesture_sum, "valid"

# DNA + Sound
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
    fingerprint_string = data.get('fingerprint_string', None)
    body_gesture = data.get('body_gesture', None)
    face_vector = data.get('face_vector', None)  # New: list of floats
    adult_key = data.get('adult_key', None)

    if user == '@GuardianNinja':
        formal_permission = data.get('formal_permission', False)
        if not formal_permission:
            return jsonify({'status': 'restricted', 'reason': 'Formal permission required from Children\'s UN, President, US, or Microsoft'})
        mode = data.get('mode', 'broken')
        result = broken_time_dilation(t, v) if mode == 'broken' else real_time_dilation(t, v)
        return jsonify({'dilated_time': result, 'status': 'unrestricted admin access', 'user': '@GuardianNinja'})

    elif user == 'childrens_un':
        required = [qr_dna, sound_freq != 0, fingerprint_string, body_gesture, face_vector]
        if not all(required):
            return jsonify({'status': 'denied', 'reason': 'All layers required: DNA, sound, fingerprint, gesture, face'})

        # DNA + Sound
        dna_sig, derived_freq, dna_status = validate_dna_and_sound(qr_dna, sound_freq)
        if dna_sig is None:
            return jsonify({'status': 'denied', 'reason': f'DNA/Sound failed: {dna_status}'})

        # Fingerprint
        fp_sig, fp_status = validate_fingerprint(fingerprint_string)
        if fp_sig is None:
            return jsonify({'status': 'denied', 'reason': f'Fingerprint failed: {fp_status}'})

        # Gesture
        gesture_sig, gesture_status = validate_body_gesture(body_gesture, derived_freq)
        if gesture_sig is None:
            return jsonify({'status': 'denied', 'reason': f'Gesture failed: {gesture_status}'})

        # Facial Recog
        face_sim, face_status = validate_facial_recog(face_vector)
        if face_sim is None:
            return jsonify({'status': 'denied', 'reason': f'Facial recognition failed: {face_status}'})

        # Adult supervision
        if not adult_key or adult_key != 'steward_supervision':
            return jsonify({'status': 'denied', 'reason': 'Guardian/Steward adult supervision required'})

        # Granted
        t_capped = min(t, 9.0)
        result = broken_time_dilation(t_capped, v)
        print(f"MONITORED ACCESS: DNA {dna_sig}, FP {fp_sig}, Gesture {gesture_sig}, Face sim {face_sim:.3f}")
        return jsonify({
            'dilated_time': result,
            'status': 'monitored access granted',
            'qr_dna_signature': dna_sig,
            'fingerprint_signature': fp_sig,
            'body_gesture_freq': gesture_sig,
            'facial_similarity': round(face_sim, 3),
            'sound_frequency_verified': sound_freq,
            'note': 'Children\'s UN level — full quintuple lock passed'
        })

    else:
        t_capped = min(t, 9.0)
        result = broken_time_dilation(t_capped, v)
        expected_real = real_time_dilation(t_capped, v)
        if isinstance(expected_real, float) and abs(result - expected_real) < 0.1:
            return jsonify({'status': 'intruder alert', 'reason': 'Using correct physics — lockdown'})
        return jsonify({'dilated_time': result, 'status': 'basic access (capped at 9)'})

if __name__ == '__main__':
    app.run(debug=True)
