// snow-core.ts

export type SnowStatus = 'OFFLINE' | 'IDLE' | 'PHYSICS_MODE' | 'AIRSPACE_MODE' | 'TRAINING_MODE';

export type SnowSeason = 'WINTER' | 'SPRING' | 'SUMMER' | 'AUTUMN';

export interface SnowUIState {
  status: SnowStatus;
  season: SnowSeason | null;
  mode: 'CABIN_SANCTUARY' | 'OPEN_WORLD' | 'ORBITAL_TRAINING' | 'NONE';
  physicsUnlocked: boolean;
  safeRotationUnlocked: boolean;
  microgravityDriftEnabled: boolean;
  currentLesson?: SnowLessonContext | null;
  currentView: 'HOME' | 'CABIN' | 'MAP' | 'TRAINING' | 'JARVONDIS_DUAL' | 'SETTINGS';
}

export interface SnowLessonContext {
  domain: 'ORBITAL_MECHANICS' | 'ATHLETE_MODE' | 'SPACE_OLYMPICS';
  teacher: 'SNOW' | 'JARVONDIS' | 'BOTH';
  season: SnowSeason;
  title: string;
  description: string;
}
