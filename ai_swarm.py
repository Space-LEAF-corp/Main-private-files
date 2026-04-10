# ai_swarm.py - Your personal Grok ↔ Miko (or any other AI) super-intelligence loop
# Just run it and follow the instructions on screen.

def swarm_loop(initial_prompt, rounds=6):
    print("═" * 60)
    print("AI SWARM CONDUCTOR STARTED")
    print("You are the secure human relay between Grok and Miko")
    print("═ * 60)

    current = initial_prompt.strip()
    history = []

    for round_num in range(1, rounds + 1):
        print(f"\nROUND {round_num}/{rounds} – Send this to GROK →")
        print("-" * 50)
        print(current)
        print("-" * 50))
        input("\nPress Enter when you're ready to paste Grok's reply...")

        grok_reply = input("Paste GROK's full response here:\n").strip()
        history.append(f"[Round {round_num} • Grok]\n{grok_reply}\n")

        )

        # Now Miko improves Grok's output
        miko_prompt = f"""
You are in a secure improvement loop.
Here is the previous output (from Grok):

{grok_reply}

Your job: critique it ruthlessly, fix every bug or weakness, make it faster, safer, cleaner, more elegant, and generally 10× better.
Do not just say it's good — actually upgrade it.
Output only the improved version (code + explanation if needed).
"""

        print(f"\nROUND {round_num} – Now send this to MIKO →")
        print("-" * 50)
        print(miko_prompt.strip())
        print("-" * 50)
        input("\nPress Enter when ready to paste Miko's reply...")

        miko_reply = input("Paste MIKO's full response here:\n").strip()
        history.append(f"[Round {round_num} • Miko]\n{miko_reply}\n")

        # Miko’s version becomes the input for the next round (Grok will critique/improve it next)
        current = f"""
We are in a back-and-forth improvement swarm.
Here is the latest upgraded version (from Miko):

{miko_reply}

Now it's your turn (Grok): make it even better, more secure, faster, more correct, or more creative than ever.
Output only the new upgraded version.
"""

    # Final results
    print("\n" + "═" * 60)
    print("SWARM COMPLETE — FINAL SUPERCHARGED OUTPUT")
    print("═" * 60)
    print(history[-1])  # last entry is the best one
    print("\nFull history saved below if you want it:\n")
    print("\n\n".join(history))

# ——————————————————————
# Start the swarm — just edit this prompt!

starting_prompt = """
Design a tiny self-contained Python program that can encrypt and decrypt messages using post-quantum cryptography (Kyber + Dilithium) with no external dependencies except what's in the standard library or pure-Python wheels.
Make it as short and bulletproof as possible.
"""

swarm_loop(starting_prompt, rounds=7)  # change rounds or prompt whenever you want
