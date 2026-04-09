import random
import time

class SafetyValidator:
    """Validates the composition for safety: ensures only food-safe, natural/synthetic elements that pass through without alteration."""
    SAFE_NATURAL_ELEMENTS = ['cellulose', 'pectin', 'chlorophyll', 'vitamins', 'minerals', 'herbal extracts']
    SAFE_SYNTHETIC_ELEMENTS = ['biodegradable polymers', 'inert silicon microsensors', 'edible dyes']
    
    def __init__(self, composition):
        self.composition = composition
    
    def is_safe(self):
        """Check if all elements are safe, excretable, and non-altering."""
        for element in self.composition:
            if element not in self.SAFE_NATURAL_ELEMENTS + self.SAFE_SYNTHETIC_ELEMENTS:
                raise ValueError(f"Unsafe element detected: {element}. Must be natural/synthetic and food-safe.")
        return True

class HeartShapedHerb:
    """Observable object tool: Simulates a safe, ingestible sensor that observes internal systems without alteration."""
    
    def __init__(self, composition, observer_callback=None):
        """Initialize with composition and optional observer (e.g., app notifier)."""
        self.validator = SafetyValidator(composition)
        self.validator.is_safe()  # Enforce safety at init
        self.composition = composition
        self.position = 'pre-ingestion'
        self.data_log = {}  # Collected observational data
        self.observer_callback = observer_callback  # For real-time updates (observable pattern)
        self.excreted = False
    
    def ingest(self):
        """Simulate ingestion: Start the pass-through process."""
        self.position = 'mouth/esophagus'
        self._notify_observer(f"Ingested. Composition validated as safe: {self.composition}")
        self._collect_data(self.position)
    
    def pass_through(self, organ):
        """Simulate transit to next organ/system, observing without alteration."""
        if self.excreted:
            raise ValueError("Herb has already been excreted.")
        self.position = organ
        self._notify_observer(f"Passing through {organ}. No alterations detected.")
        time.sleep(0.5)  # Simulate transit time
        self._collect_data(organ)
    
    def _collect_data(self, location):
        """Simulate microscopic/atomic-level scan: Passive observation of environmental/anatomical data."""
        # Randomized but realistic data for simulation (e.g., pH, temp, structural metrics)
        ph = round(random.uniform(1.5, 7.5), 2)  # e.g., stomach acid vs. intestines
        temp = round(random.uniform(36.5, 37.5), 1)  # Body temp range
        moisture = round(random.uniform(70, 95), 0)  # % humidity in tract
        # Pseudo-atomic structure: Simulate counts of atoms/molecules (non-invasive, observational only)
        atomic_structure = {
            'carbon_atoms': random.randint(1000, 5000),
            'hydrogen_atoms': random.randint(2000, 10000),
            'oxygen_atoms': random.randint(500, 3000),
            'cellular_density': round(random.uniform(1e6, 1e9), 0),  # Cells per mm3
            'microbiome_signature': random.choice(['balanced', 'acidic', 'alkaline']),  # Environmental flag
        }
        self.data_log[location] = {
            'pH': ph,
            'temperature_C': temp,
            'moisture_%': moisture,
            'atomic_level_structure': atomic_structure,
            'flow_rate': round(random.uniform(0.1, 1.0), 2),  # cm/s simulation
            'texture': random.choice(['smooth', 'mucosal', 'contractile']),
            'notes': 'Observation only: No interaction or alteration.'
        }
        self._notify_observer(f"Data collected at {location}: {self.data_log[location]}")
    
    def excrete(self):
        """Simulate natural excretion: End process and return full data for professional review."""
        self.position = 'excreted'
        self.excreted = True
        self._notify_observer("Excreted safely. Full data catalog ready for interpretation.")
        return self.data_log
    
    def _notify_observer(self, message):
        """Observable pattern: Notify callback (e.g., app) with updates."""
        if self.observer_callback:
            self.observer_callback(message)

# Example usage: Simulate the full process with an app notifier​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​​
