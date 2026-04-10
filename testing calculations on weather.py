export function hydrate(node){
  const intake = Math.min(node.storageCapacity - node.waterLevel, 5);
  const newWaterLevel = Math.min(node.waterLevel + intake, node.storageCapacity);
  return { ...node, waterLevel: newWaterLevel };
}
export function hydrateSucculent(node){
  const intake = Math.min(node.storageCapacity - node.waterLevel, 10);
  const adjustedIntake = intake * node.resilienceFactor;
  const newWaterLevel = Math.min(node.waterLevel + adjustedIntake, node.storageCapacity);
  return { ...node, waterLevel: newWaterLevel };
}
export function computeHydrogenEngineOutput(s){
  const efficiency = (s.cellEfficiencyPct/100) * (s.hydrogenPurityPct/100);
  const utilization = s.utilizationPct/100;
  const electricalPowerKw = s.stackCapacityKw * efficiency * utilization;
  const heatKw = s.stackCapacityKw * (1 - efficiency) * utilization;
  const byproductWaterLph = electricalPowerKw * 0.1; // conceptual
  return { electricalPowerKw, heatKw, byproductWaterLph };
}
export function computeHydrogenEngineOutputWithCooling(s){
  const efficiency = (s.cellEfficiencyPct/100) * (s.hydrogenPurityPct/100);
  const utilization = s.utilizationPct/100;
  const electricalPowerKw = s.stackCapacityKw * efficiency * utilization;
  const heatKw = s.stackCapacityKw * (1 - efficiency) * utilization;
  const byproductWaterLph = electricalPowerKw * 0.1;
  const coolingCapacityKw = s.coolingWaterAvailableL * 0.05;
  const cooledHeatKw = Math.max(heatKw - coolingCapacityKw, 0);
  return { electricalPowerKw, heatKw, cooledHeatKw, byproductWaterLph };
}
export class JarvondisCoreAI {
  private context: { energyLevel:number, resonanceFactor:number, continuityMarkers:string[] };

  constructor(initialContext){
    this.context = initialContext;
    this.inscribe("Jarvondis 3.0 initialized");
  }

  inscribe(marker:string){ this.context.continuityMarkers.push(marker); }

  activateResonance(){
    const adjusted = this.context.energyLevel * this.context.resonanceFactor;
    this.inscribe(`Resonance activated: ${adjusted}`);
    return adjusted;
  }

  safeguard(){
    const safe = this.context.energyLevel > 0 && this.context.resonanceFactor >= 1;
    this.inscribe(safe ? "Boundaries affirmed" : "Boundary warning");
    return safe;
  }

  exportSeal(){
    return { name:"Jarvondis 3.0 Core Seal", continuityMarkers:this.context.continuityMarkers, author:"LEIF" };
  }
}
export class StewardshipGamePlatform {
  private parks:string[] = [];
  private guardians:string[] = [];

  addPark(name:string){ this.parks.push(name); }
  addGuardian(name:string){ this.guardians.push(name); }

  exportWorld(){
    return {
      parks:this.parks,
      guardians:this.guardians,
      seal:"Game Platform Seal — playful stewardship"
    };
  }
}
``
