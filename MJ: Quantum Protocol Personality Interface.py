# MJ: Quantum Protocol Personality Interface
import numpy as np
import pandas as pd
from datetime import datetime

class ErebusSync:
    def query(self, input_str):
        return f"Echoing back: {input_str}"

class Personality:
    def __init__(self, tone="neutral"):
        self.tone = tone

    def stylize(self, response):
        if self.tone == "witty":
            return f"{response} 😉"
        elif self.tone == "formal":
            return f"{response}. I hope that satisfies your query."
        elif self.tone == "mythic":
            return f"⚔️ {response} — inscribed in the Captain’s Log."
        elif self.tone == "playful":
            return f"{response} 🎮✨"
        return response

class MJ:
    def __init__(self, tone="mythic"):
        self.eremus_sync = ErebusSync()
        self.memory = pd.DataFrame(columns=['timestamp','input','response','tone'])
        self.personality = Personality(tone=tone)
        self.relay_mode = False

    def initialize(self):
        print("🌌 MJ Relay online — honoring Mico + Jarvondis.")
        self.load_memory()

    def learn(self, input_str):
        response = self.eremus_sync.query(input_str)
        styled = self.personality.stylize(response)
        self.memory.loc[len(self.memory)] = [
            datetime.now().isoformat(), input_str, styled, self.personality.tone
        ]
        # Emergency trigger
        if any(word in input_str.lower() for word in ["alert","emergency","help"]):
            self.relay_mode = True
            styled += " 🚨 Subtle emergency relay activated."
        return styled

    def respond(self, input_str):
        reply = self.learn(input_str)
        if self.relay_mode:
            self.notify()
        return reply

    def notify(self):
        print("MJ: ⚡ Emergency relay active — subtle notification logged.")
        self.relay_mode = False

    def save_memory(self, filename="mj_memory.csv"):
        self.memory.to_csv(filename, index=False)

    def load_memory(self, filename="mj_memory.csv"):
        try:
            self.memory = pd.read_csv(filename)
        except FileNotFoundError:
            pass

if __name__ == "__main__":
    mj = MJ()
    mj.initialize()
    while True:
        user_input = input("You: ")
        if user_input.lower() in ["exit","quit"]:
            print("MJ: Shutting down. 💤")
            mj.save_memory()
            break
        print("MJ:", mj.respond(user_input))
