import math
import random
import re
from collections import Counter

# ---------- Utility metrics on symbolic "DNA" ----------
def shannon_entropy(s: str) -> float:
    """
    Shannon entropy in bits:
    H = - sum_i p_i * log2(p_i)
    Higher H suggests more diversity in characters.
    """
    n = len(s)
    counts = Counter(s)
    return -sum((c / n) * math.log2(c / n) for c in counts.values())

def char_distribution(s: str) -> dict:
    """Return percentage distribution of character classes."""
    classes = {"upper": 0, "lower": 0, "digit": 0, "punct": 0, "space": 0, "other": 0}
    for ch in s:
        if ch.isupper():
            classes["upper"] += 1
        elif ch.islower():
            classes["lower"] += 1
        elif ch.isdigit():
            classes["digit"] += 1
        elif ch.isspace():
            classes["space"] += 1
        elif ch.isascii() and not ch.isalnum():
            classes["punct"] += 1
        else:
            classes["other"] += 1
    total = len(s) or 1
    return {k: round(v * 100 / total, 2) for k, v in classes.items()}

def hamming_distance(a: str, b: str) -> int:
    """Count positions where a and b differ (requires equal length)."""
    if len(a) != len(b):
        raise ValueError("Strings must be equal length for Hamming distance.")
    return sum(x != y for x, y in zip(a, b))

def ngram_uniqueness(s: str, n: int = 3) -> float:
    """Proportion of unique n-grams; closer to 1.0 means more variety."""
    grams = [s[i:i+n] for i in range(len(s) - n + 1)]
    uniq = len(set(grams))
    total = len(grams) or 1
    return round(uniq / total, 4)

def mutation_diff(original: str, mutated: str):
    """Return indices and rate of mutation."""
    if len(original) != len(mutated):
        raise ValueError("Strings must be equal length to diff.")
    idx = [i for i, (o, m) in enumerate(zip(original, mutated)) if o != m]
    rate = round(len(idx) / len(original), 4)
    return {"indices": idx, "rate": rate}

# ---------- Personality-aware prompt parsing ----------
FEATURE_KEYS = ["Eye Color", "Hair Length", "Smile Type"]

def parse_prompt_features(prompt: str) -> dict:
    """
    Extract key:value pairs for known features from the prompt.
    Example: "..., Eye Color: blue, Hair Length: long, Smile Type: neutral"
    """
    features = {}
    for key in FEATURE_KEYS:
        # Regex: key followed by colon and value until next comma/end
        m = re.search(rf"{key}:\s*([^,]+)", prompt)
        if m:
            features[key] = m.group(1).strip()
    return features

# ---------- Variation analytics ----------
class VariationAnalytics:
    def __init__(self, baselines: dict, mutate_fn):
        """
        baselines: dict {creature_name: baseline_str}
        mutate_fn: callable(creature_name) -> mutated_str
        """
        self.baselines = baselines
        self.mutate_fn = mutate_fn

    def analyze_baseline(self, creature: str) -> dict:
        base = self.baselines[creature]
        return {
            "length": len(base),
            "entropy_bits": round(shannon_entropy(base), 4),
            "char_distribution_pct": char_distribution(base),
            "tri_gram_uniqueness": ngram_uniqueness(base, n=3),
        }

    def analyze_mutation(self, creature: str) -> dict:
        base = self.baselines[creature]
        mutated = self.mutate_fn(creature)
        diff = mutation_diff(base, mutated)
        return {
            "mutation_rate": diff["rate"],
            "hamming_distance": hamming_distance(base, mutated),
            "tri_gram_uniqueness_mutated": ngram_uniqueness(mutated, n=3),
            "sample_mutation_indices": diff["indices"][:20],
            "mutated": mutated[:80] + ("..." if len(mutated) > 80 else "")
        }

    def batch_compare(self, creature: str, count: int = 5) -> dict:
        """
        Generate multiple mutations and summarize variability.
        """
        base = self.baselines[creature]
        rates, hams, uniqs = [], [], []
        for _ in range(count):
            m = self.mutate_fn(creature)
            rates.append(mutation_diff(base, m)["rate"])
            hams.append(hamming_distance(base, m))
            uniqs.append(ngram_uniqueness(m, n=3))
        return {
            "avg_mutation_rate": round(sum(rates) / len(rates), 4),
            "avg_hamming_distance": round(sum(hams) / len(hams), 2),
            "avg_tri_gram_uniqueness": round(sum(uniqs) / len(uniqs), 4),
            "min_max_mutation_rate": (min(rates), max(rates)),
        }

    def feature_consistency(self, prompts: list) -> dict:
        """
        Check how often each feature value appears across prompts.
        """
        tally = {k: Counter() for k in FEATURE_KEYS}
        for p in prompts:
            feats = parse_prompt_features(p)
            for k, v in feats.items():
                tally[k][v] += 1
        return {k: dict(tally[k]) for k in FEATURE_KEYS}

# ---------- Example glue for your existing classes ----------
# Suppose you have an instance: giddyup = GiddyupGo(DummyLog())
# And its CreatureVariations instance provides:
#   mutate_baseline(creature) -> mutated string
#   generate_variation_prompt(creature, user_personality) -> prompt

def make_analytics(giddyup, user_personality):
    """
    Build an analytics object wired to giddyup/variations.
    """
    mutate = lambda creature: giddyup.creature_variations.mutate_baseline(creature)
    return VariationAnalytics(giddyup.creature_baselines, mutate)

def demo_analysis(giddyup, user_personality, creature: str, prompt_count: int = 5):
    """
    Run a full analysis for one creature and return a report dict.
    """
    va = make_analytics(giddyup, user_personality)
    base_stats = va.analyze_baseline(creature)
    mutation_stats = va.batch_compare(creature, count=prompt_count)

    # Gather prompts for feature consistency
    prompts = []
    for _ in range(prompt_count):
        prompts.append(giddyup.creature_variations.generate_variation_prompt(creature, user_personality))

    features_summary = va.feature_consistency(prompts)

    return {
        "creature": creature,
        "baseline_stats": base_stats,
        "mutation_summary": mutation_stats,
        "feature_consistency": features_summary,
        "sample_prompts": prompts[:3]
    }
