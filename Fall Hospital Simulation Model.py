# Fall Hospital Simulation Model
# Purpose: Provide a scaffold for simulating hospital operations during seasonal shifts.
# Note: This is a symbolic model, not tied to real patient data.

import random
from dataclasses import dataclass, field
from typing import List, Dict

@dataclass
class Patient:
    id: int
    condition: str
    severity: int  # 1 = mild, 5 = critical
    admitted: bool = False

@dataclass
class Staff:
    id: int
    role: str  # "doctor", "nurse", "support"
    available: bool = True

@dataclass
class Hospital:
    name: str
    capacity: int
    patients: List[Patient] = field(default_factory=list)
    staff: List[Staff] = field(default_factory=list)
    logs: List[str] = field(default_factory=list)

    def admit_patient(self, patient: Patient):
        if len(self.patients) < self.capacity:
            patient.admitted = True
            self.patients.append(patient)
            self.logs.append(f"Patient {patient.id} admitted with {patient.condition}")
        else:
            self.logs.append(f"Patient {patient.id} could not be admitted (capacity full)")

    def assign_staff(self):
        for patient in self.patients:
            available_staff = [s for s in self.staff if s.available]
            if available_staff:
                staff_member = random.choice(available_staff)
                staff_member.available = False
                self.logs.append(f"Staff {staff_member.id} ({staff_member.role}) assigned to Patient {patient.id}")
            else:
                self.logs.append(f"No staff available for Patient {patient.id}")

    def discharge_patient(self, patient_id: int):
        for patient in self.patients:
            if patient.id == patient_id:
                self.patients.remove(patient)
                self.logs.append(f"Patient {patient.id} discharged")
                break

    def simulate_day(self):
        # Random admissions
        new_patients = [Patient(id=random.randint(1000, 9999),
                                condition=random.choice(["flu", "injury", "respiratory"]),
                                severity=random.randint(1, 5)) for _ in range(random.randint(1, 5))]
        for p in new_patients:
            self.admit_patient(p)

        # Assign staff
        self.assign_staff()

        # Random discharges
        if self.patients:
            discharged = random.choice(self.patients)
            self.discharge_patient(discharged.id)

        self.logs.append("Day simulation complete")

# Example usage
hospital = Hospital(name="Fall General", capacity=10,
                    staff=[Staff(id=i, role=random.choice(["doctor", "nurse", "support"])) for i in range(1, 6)])

for day in range(3):  # simulate 3 days
    hospital.simulate_day()

print("Simulation Logs:")
for log in hospital.logs:
    print(log)
