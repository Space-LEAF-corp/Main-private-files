// Carbon puck core types

export type PuckState =
  | "IN_VEHICLE"
  | "AT_SERVICE_CENTER"
  | "IN_TRANSIT"
  | "AT_LAB"
  | "PROCESSED"
  | "DISPOSED_OR_STORED";

export interface CustodyEvent {
  eventId: string;
  puckId: string;
  timestamp: string; // ISO 8601
  locationCode: string; // e.g. ISO country/region or facility code
  actorRealId: string; // hashed or pseudonymous Real ID
  actorRole: "service_tech" | "transporter" | "lab_tech" | "regulator";
  fromState: PuckState;
  toState: PuckState;
  notes?: string;
  eventSignature: string; // cryptographic signature
}

export interface CarbonPuck {
  puckId: string;
  vehicleIdHash: string;
  ownerIdHash?: string;

  massGrams: number;
  estimatedCarbonFraction: number; // 0–1

  estimatedCo2EquivalentGrams: number;

  engineType: "diesel" | "gasoline" | "hybrid" | "other";
  fuelType: string;
  regionCode: string;

  captureTimestamp: string;
  removalTimestamp: string;

  createdByDeviceId: string;
  createdByRealId: string;
  digitalSignature: string;

  custodyEvents: CustodyEvent[];
}

// Simple CO2-equivalent calculation
export function computeCo2EquivalentGrams(
  massGrams: number,
  estimatedCarbonFraction: number
): number {
  const carbonMass = massGrams * estimatedCarbonFraction;
  const co2Mass = carbonMass * (44 / 12);
  return co2Mass;
}

// Example Leaf Coin issuance logic
export interface LeafIssuanceResult {
  puckId: string;
  co2EquivalentGrams: number;
  leafCoinsIssued: number;
  reason: string;
}

export function issueLeafCoinsForPuck(
  puck: CarbonPuck,
  gramsPerLeaf: number = 1000
): LeafIssuanceResult {
  const co2 = puck.estimatedCo2EquivalentGrams;

  // Example rule: only issue if puck reached a terminal state
  const lastEvent = puck.custodyEvents[puck.custodyEvents.length - 1];
  const terminalStates: PuckState[] = ["PROCESSED", "DISPOSED_OR_STORED"];

  if (!lastEvent || !terminalStates.includes(lastEvent.toState)) {
    return {
      puckId: puck.puckId,
      co2EquivalentGrams: co2,
      leafCoinsIssued: 0,
      reason: "Puck has not reached a terminal custody state."
    };
  }

  const coins = Math.floor(co2 / gramsPerLeaf);

  return {
    puckId: puck.puckId,
    co2EquivalentGrams: co2,
    leafCoinsIssued: coins,
    reason: "Puck processed with valid terminal state."
  };
}
