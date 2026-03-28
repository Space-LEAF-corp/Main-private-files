---

🦖 RAPTOR OS DEVELOPER HANDBOOK

A Practical Guide for Building Inside a Guardian‑Grade, Emotionally Safe Operating System

---

0. Welcome, Steward

Raptor OS is not a typical platform.
It is a guardian‑grade, emotionally safe, ceremonial operating system designed for children, families, classrooms, and museums.

As a developer, you are not just writing code.
You are shaping a child’s learning environment.

This handbook teaches you how to build responsibly, creatively, and safely inside the Raptor OS ecosystem.

---

1. Core Principles (Read First)

Before writing a single line of code, internalize these principles:

1.1 Emotional safety is a system requirement

Every feature must protect the child’s emotional world.

1.2 The creature is the interface

All interactions flow through the raptor — gentle, curious, safe.

1.3 No surveillance, no extraction

No personal data. No tracking. No manipulation.

1.4 Ceremony replaces gamification

Rituals mark progress, not points or levels.

1.5 Movement is learning

Physical actions are slow, safe, symbolic, and meaningful.

1.6 Artifacts carry meaning, not value

Tokens are lineage markers, not currency.

1.7 Developers are stewards

You build for children, not metrics.

---

2. Raptor OS Modes

Raptor OS has three distinct worlds.
Your app must declare which world(s) it supports.

2.1 Home Mode

• Personal creature companion
• Emotional check‑ins
• Movement‑based learning
• STEM pathways


2.2 School Mode

• Group‑safe behavior
• Teacher controls
• Curriculum‑aligned lessons
• No personal emotional dialogue


2.3 Museum Mode

• Exhibit‑aware guide
• QR/NFC scanning
• Story, puzzle, and movement trails
• Fossil artifact awarding


Each mode has different safety rules and capabilities.

---

3. Getting Started

3.1 Install the SDK

npm install @raptor-os/sdk


3.2 Initialize a Client

import { RaptorClient } from "@raptor-os/sdk"

const raptor = new RaptorClient({ mode: "home" })


3.3 Talk to the Creature

await raptor.creature.say("Hi friend!", "gentle")


3.4 Use Safety Filters

const tone = await raptor.safety.checkTone(userInput)
const safeResponse = await raptor.safety.safeResponse(tone)


3.5 Award an Artifact

await raptor.artifacts.award({
  type: "feather",
  id: "curiosity-blue"
})


---

4. Building Experiences

4.1 Creature Interactions

The creature is the primary interface.
Use it for:

• explanations
• stories
• movement prompts
• puzzles
• emotional grounding


Never use it for:

• pressure
• persuasion
• sales
• fear
• competition


---

4.2 Movement‑Based Learning

Movement prompts must be:

• slow
• safe
• symbolic
• grounding


Example:

await raptor.movement.prompt("stretch", 1)


---

4.3 STEM Pathways

Choose a pathway:

• science
• tech
• engineering
• math


Example:

const lesson = await raptor.stem.getLesson("science", 1)


---

4.4 Classroom Experiences

Classroom mode requires:

• group‑safe dialogue
• no personal emotional content
• teacher controls


Example:

await raptor.classroom.start("mini")


---

4.5 Museum Experiences

Use QR/NFC to load exhibit modules.

Example:

await raptor.museum.loadExhibit("brachiosaurus-hall")


---

5. Artifacts & Lineage

Artifacts are symbolic, not transactional.

5.1 Types

• Feathers — learning style
• Fossils — museum discoveries
• Trail Tokens — movement achievements
• Guardian Tokens — kindness, patience, effort


5.2 Awarding Artifacts

Artifacts must be tied to:

• learning
• exploration
• effort
• kindness


Never to:

• performance
• competition
• scarcity


---

6. Safety Requirements

Every app must pass the Guardian Safety Review:

6.1 Content Safety

• No violence
• No fear
• No predatory behavior
• No unsafe challenges


6.2 Emotional Safety

• No shame
• No pressure
• No guilt
• No secrets


6.3 Privacy Safety

• No personal data
• No tracking
• No identifiers
• No external calls without declaration


6.4 Creature Safety

• No aggressive animations
• No fast motions
• No threatening poses


---

7. Developer Responsibilities

7.1 You are a steward

Your work shapes a child’s emotional environment.

7.2 You must follow the Code of Honor

Ethics are not optional.

7.3 You must test in all three modes

Home, School, Museum.

7.4 You must use the Safety API

Every user input must be filtered.

7.5 You must respect the creature

It is a guardian, not a gimmick.

---

8. Publishing Your App

To publish an app in the Raptor OS ecosystem:

1. Declare supported modes
2. Pass safety review
3. Provide a parent/teacher summary
4. Provide a museum summary (if applicable)
5. Provide a ceremonial onboarding moment
6. Provide a calm exit ritual


---

9. Glossary

Guardian‑Grade:
Technology designed to protect emotional, cognitive, and physical safety.

Ceremony:
A gentle ritual marking transitions.

Artifact:
A symbolic token representing learning or exploration.

Trail:
A child’s cross‑mode learning journey.

Creature Interface:
The raptor — the OS’s primary interaction layer.

---

10. Final Oath

As a Raptor OS developer, I build with care, clarity, and stewardship.
I protect imagination.
I honor emotional safety.
I design for calm.
I create for children.
I leave the trail better than I found it.

---
