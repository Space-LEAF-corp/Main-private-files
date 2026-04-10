// models.ts
export type AgeBand = "5-7" | "8-10" | "11-13";

export type FactCategory =
  | "Nature"
  | "Science"
  | "Math"
  | "History"
  | "Arts"
  | "Feelings"
  | "Community"
  | "HealthyHabits";

export interface FactCreature {
  id: string;                 // e.g., fc_nature_ants_001
  name: string;               // e.g., "Ant Scout"
  category: FactCategory;
  ageBand: AgeBand;
  prompt: string;             // child-facing discovery prompt
  evidenceType: "Text" | "Drawing" | "Photo" | "Count" | "Audio"; // choose what’s allowed
  captureRule: {              // what counts as a valid capture
    minWords?: number;
    minCount?: number;
    requiresObservation?: boolean;
  };
  parentNote?: string;        // auto-generated teaching note for parents
}

export interface CaptureEvidence {
  type: "Text" | "Drawing" | "Photo" | "Count" | "Audio";
  text?: string;
  count?: number;
  drawingId?: string;
  photoId?: string;
  audioId?: string;
}

export interface Capture {
  captureId: string;
  creatureId: string;
  childId: string;            // pseudonymous local ID
  timestamp: string;          // ISO string
  evidence: CaptureEvidence;
  verified: boolean;          // local “grownup check” toggle
}

export interface DailyJournal {
  childId: string;
  date: string;
  captures: Capture[];
  feelingsNote?: string;      // optional reflection
}
[
  {
    "id": "fc_nature_clouds_001",
    "name": "Cloud Whisper",
    "category": "Nature",
    "ageBand": "5-7",
    "prompt": "Find a cloud. What shape does it look like?",
    "evidenceType": "Text",
    "captureRule": { "minWords": 3, "requiresObservation": true },
    "parentNote": "Child practiced observation and metaphor using weather as a prompt."
  },
  {
    "id": "fc_math_count_rocks_001",
    "name": "Rock Counter",
    "category": "Math",
    "ageBand": "5-7",
    "prompt": "Count small rocks you see today. How many did you find?",
    "evidenceType": "Count",
    "captureRule": { "minCount": 1, "requiresObservation": true },
    "parentNote": "Counting practice and one-to-one correspondence in a real context."
  },
  {
    "id": "fc_science_ants_001",
    "name": "Ant Scout",
    "category": "Science",
    "ageBand": "8-10",
    "prompt": "Spot ants at work. What are they carrying or doing?",
    "evidenceType": "Text",
    "captureRule": { "minWords": 5, "requiresObservation": true },
    "parentNote": "Introduces social insects and task specialization."
  },
  {
    "id": "fc_community_kindness_001",
    "name": "Kindness Lantern",
    "category": "Community",
    "ageBand": "5-7",
    "prompt": "Did you help someone today? Describe your kind action.",
    "evidenceType": "Text",
    "captureRule": { "minWords": 5 },
    "parentNote": "Encourages empathy and social reflection."
  }
]
// capture.ts
import { Capture, CaptureEvidence, FactCreature } from "./models";

export function validateEvidence(creature: FactCreature, evidence: CaptureEvidence): string[] {
  const errors: string[] = [];
  if (evidence.type !== creature.evidenceType) {
    errors.push("Evidence type does not match this creature’s prompt.");
  }
  if (creature.captureRule.minWords && (!evidence.text || wordCount(evidence.text) < creature.captureRule.minWords)) {
    errors.push(`Please write at least ${creature.captureRule.minWords} words.`);
  }
  if (creature.captureRule.minCount && (!evidence.count || evidence.count < creature.captureRule.minCount)) {
    errors.push(`Please enter a count of at least ${creature.captureRule.minCount}.`);
  }
  return errors;
}

function wordCount(s: string) {
  return s.trim().split(/\s+/).filter(Boolean).length;
}

export function createCapture(
  creature: FactCreature,
  childId: string,
  evidence: CaptureEvidence
): { capture?: Capture; errors?: string[] } {
  const errors = validateEvidence(creature, evidence);
  if (errors.length) return { errors };

  const capture: Capture = {
    captureId: `cap_${Date.now()}_${Math.random().toString(36).slice(2)}`,
    creatureId: creature.id,
    childId,
    timestamp: new Date().toISOString(),
    evidence,
    verified: false
  };
  return { capture };
}

export function verifyCapture(capture: Capture, grownupApproved: boolean): Capture {
  return { ...capture, verified: grownupApproved };
}
// summary.ts
import { DailyJournal, FactCreature, Capture } from "./models";

export function buildParentSummary(
  journal: DailyJournal,
  creatureIndex: Map<string, FactCreature>
) {
  const items = journal.captures.map((cap: Capture) => {
    const creature = creatureIndex.get(cap.creatureId);
    return {
      name: creature?.name ?? "Unknown",
      category: creature?.category ?? "Unknown",
      prompt: creature?.prompt ?? "",
      childEvidence: cap.evidence,
      verifiedByGrownup: cap.verified,
      parentTeachingNote: creature?.parentNote ?? ""
    };
  });

  return {
    date: journal.date,
    feelingsNote: journal.feelingsNote ?? "",
    items
  };
}
// CaptureCard.tsx
import React, { useState } from "react";
import { FactCreature, CaptureEvidence } from "./models";
import { createCapture } from "./capture";

export const CaptureCard: React.FC<{
  creature: FactCreature;
  childId: string;
  onCaptured: (cap: any) => void;
}> = ({ creature, childId, onCaptured }) => {
  const [text, setText] = useState("");
  const [count, setCount] = useState<number | undefined>();
  const [errors, setErrors] = useState<string[]>([]);

  function submit() {
    const evidence: CaptureEvidence =
      creature.evidenceType === "Text"
        ? { type: "Text", text }
        : { type: "Count", count: count ?? 0 };

    const result = createCapture(creature, childId, evidence);
    if (result.errors) setErrors(result.errors);
    else if (result.capture) onCaptured(result.capture);
  }

  return (
    <div className="card">
      <h3>{creature.name}</h3>
      <p>{creature.prompt}</p>

      {creature.evidenceType === "Text" && (
        <textarea value={text} onChange={(e) => setText(e.target.value)} placeholder="Write what you observed..." />
      )}

      {creature.evidenceType === "Count" && (
        <input type="number" value={count ?? ""} onChange={(e) => setCount(Number(e.target.value))} placeholder="Enter a count" />
      )}

      {errors.length > 0 && <ul>{errors.map((e) => <li key={e}>{e}</li>)}</ul>}

      <button onClick={submit}>Capture</button>
    </div>
  );
};
