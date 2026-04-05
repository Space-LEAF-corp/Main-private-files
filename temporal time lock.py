import random

class Whiskers:
    def say(self, text):
        print(f"[Whiskers] {text}")

whiskers = Whiskers()

# -----------------------------
# Dice utility
# -----------------------------
def roll_die():
    return random.randint(1, 6)

def get_user_dice_input(count):
    values = []
    for i in range(count):
        while True:
            try:
                val = int(input(f"Enter die {i+1} (1–6): "))
                if 1 <= val <= 6:
                    values.append(val)
                    break
                else:
                    print("Value must be between 1 and 6.")
            except ValueError:
                print("Please enter a number.")
    return values

# -----------------------------
# Standard User Flow
# -----------------------------
def standard_user_flow():
    whiskers.say("Standard user selected. Three single-die checks.")

    answers = []
    prompts = [
        "Roll one die and tell me the number.",
        "Roll again. What number do you see now.",
        "Final roll. What did you get."
    ]

    for prompt in prompts:
        whiskers.say(prompt)
        answers.append(get_user_dice_input(1)[0])

    whiskers.say(f"Standard user authentication complete. Values: {answers}")
    return answers

# -----------------------------
# Captain User Flow
# -----------------------------
def captain_user_flow():
    whiskers.say("Captain user selected. Nine questions, four requiring two dice.")

    questions = [
        ("Launch: origin marker", 1),
        ("Launch: origin lock", 2),
        ("Trajectory: path anchor", 1),
        ("Trajectory: mid-course vector", 2),
        ("Return: return marker", 1),
        ("Return: re-entry lock", 2),
        ("Final check: temporal checksum", 1),
        ("Final check: signature pair", 2),
        ("Final check: closing marker", 1)
    ]

    results = []

    for prompt, count in questions:
        whiskers.say(prompt)
        values = get_user_dice_input(count)
        results.append(values)

    whiskers.say("Captain authentication complete.")
    whiskers.say(f"Collected sequence: {results}")
    return results

# -----------------------------
# Main Menu
# -----------------------------
def main():
    whiskers.say("Welcome to the Temporal Time Lock prototype.")
    whiskers.say("Choose your mode:")
    print("1. Standard User")
    print("2. Captain User")

    while True:
        choice = input("Enter 1 or 2: ").strip()
        if choice == "1":
            standard_user_flow()
            break
        elif choice == "2":
            captain_user_flow()
            break
        else:
            print("Invalid choice. Enter 1 or 2.")

if __name__ == "__main__":
    main()