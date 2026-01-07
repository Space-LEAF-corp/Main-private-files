# Heritage Guardian Network (HGN)
# A protective, non-weaponized emergency protocol for safeguarding people
# and historic sites. Inspired by the idea of "animated guardians," but fully original.

from typing import List, Dict

class GuardianUnit:
    def __init__(self, unit_count: int = 10):
        self.unit_count = unit_count
        self.status = "dormant"
        self.units: List[str] = [f"Heritage Guardian {i+1}" for i in range(unit_count)]
    
    def activate(self, crisis_level: str) -> str:
        if self.status == "dormant":
            self.status = "active"
            return f"Activating {self.unit_count} Heritage Guardians for {crisis_level} protection."
        return "Guardians already active and deployed."
    
    def deploy(self, location: str) -> Dict[str, str]:
        if self.status == "active":
            return {
                unit: f"Deployed to {location} — providing shielding and civilian guidance."
                for unit in self.units
            }
        return {"Error": "Guardians not active yet. Activate first."}
    
    def standby(self) -> str:
        self.status = "dormant"
        return "Guardians returning to standby. Area stabilized."

class HeritageGuardianNetwork:
    def __init__(self, perimeter_km: float = 50, unit_count: int = 10):
        self.perimeter_km = perimeter_km
        self.guardians = GuardianUnit(unit_count)
        self.crisis_active = False
    
    def detect_crisis(self, description: str) -> str:
        self.crisis_active = True
        return f"Crisis detected: {description}. Initiating Heritage Guardian Protocol."
    
    def activate_protocol(self, user_id: str, location: str, crisis_description: str = None) -> Dict[str, any]:
        protocol_report = {
            "User": user_id,
            "Location": location,
            "Perimeter": f"{self.perimeter_km} km safety radius",
            "Status": "Base protocol active"
        }
        
        if crisis_description:
            crisis_status = self.detect_crisis(crisis_description)
            activation_status = self.guardians.activate(crisis_description)
            deployment = self.guardians.deploy(location)
            
            protocol_report.update({
                "CrisisStatus": crisis_status,
                "GuardianActivation": activation_status,
                "Deployments": deployment,
                "Protocol": "Heritage Guardian Network Active"
            })
        
        return protocol_report
    
    def deactivate_protocol(self) -> str:
        self.crisis_active = False
        return self.guardians.standby() + " All systems stable."

# Example usage:
if __name__ == "__main__":
    hgn = HeritageGuardianNetwork(perimeter_km=50, unit_count=20)
    
    result = hgn.activate_protocol(
        user_id="SpaceLeafCorpSteward",
        location="UK Heritage Sites",
        crisis_description="Public safety disturbance"
    )
    print(result)
    
    print(hgn.deactivate_protocol())
