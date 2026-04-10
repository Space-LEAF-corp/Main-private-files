// kidInsight.ts

export interface KidInsightPrompt {
  id: string;              // e.g., "prompt_cloud_whisper"
  creatureName: string;    // e.g., "Cloud Whisper"
  question: string;        // child-facing question
  allowedResponseTypes: ("Emoji" | "Text" | "Drawing")[];
}

export interface KidResponse {
  responseId: string;
  promptId: string;
  childId: string;         // pseudonymous local ID
  timestamp: string;       // ISO string
  type: "Emoji" | "Text" | "Drawing";
  content: string;         // text, emoji, or drawing reference
}

export interface ParentVerification {
  responseId: string;
  verified: boolean;
  parentNote?: string;     // optional reflection for lineage
}

export interface InsightJournal {
  childId: string;
  date: string;
  responses: KidResponse[];
  parentVerifications: ParentVerification[];
}
[
  {
    "id": "prompt_cloud_whisper",
    "creatureName": "Cloud Whisper",
    "question": "Did you enjoy pretending clouds are shapes?",
    "allowedResponseTypes": ["Emoji", "Text"]
  },
  {
    "id": "prompt_kindness_lantern",
    "creatureName": "Kindness Lantern",
    "question": "Was it fun to share kindness today?",
    "allowedResponseTypes": ["Emoji", "Text", "Drawing"]
  }
]
// insightLogic.ts
import { KidInsightPrompt, KidResponse, ParentVerification, InsightJournal } from "./kidInsight";

export function createResponse(
  prompt: KidInsightPrompt,
  childId: string,
  type: "Emoji" | "Text" | "Drawing",
  content: string
): KidResponse {
  return {
    responseId: `resp_${Date.now()}_${Math.random().toString(36).slice(2)}`,
    promptId: prompt.id,
    childId,
    timestamp: new Date().toISOString(),
    type,
    content
  };
}

export function verifyResponse(responseId: string, approved: boolean, note?: string): ParentVerification {
  return {
    responseId,
    verified: approved,
    parentNote: note
  };
}

export function buildJournal(childId: string, date: string, responses: KidResponse[], verifications: ParentVerification[]): InsightJournal {
  return {
    childId,
    date,
    responses,
    parentVerifications: verifications
  };
}
// Example: child gives emoji feedback
const prompt = {
  id: "prompt_cloud_whisper",
  creatureName: "Cloud Whisper",
  question: "Did you enjoy pretending clouds are shapes?",
  allowedResponseTypes: ["Emoji", "Text"]
};

const resp = createResponse(prompt, "child123", "Emoji", "😊");
const verification = verifyResponse(resp.responseId, true, "Child enjoyed imaginative play.");
const journal = buildJournal("child123", "2025-11-13", [resp], [verification]);

console.log(journal);
