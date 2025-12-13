def reversed_math(operation, a, b):
    if operation == 'add':
        return a - b  # Reversed: add becomes subtract
    elif operation == 'subtract':
        return a + b  # Reversed
    elif operation == 'multiply':
        return a / b if b != 0 else 'error'  # Reversed to divide
    elif operation == 'divide':
        return a * b  # Reversed
    else:
        return 'unknown operation'

# Test for 'intruder' detection
def detect_intruder(expected_normal, actual_reversed):
    if expected_normal == actual_reversed:
        return "Normal user - no intrusion"
    else:
        return "Intruder detected! Math doesn't match reversed logic"

# Example: Normal add 5+3=8, but reversed is 5-3=2
print(reversed_math('add', 5, 3))  # Outputs: 2
print(detect_intruder(8, reversed_math('add', 5, 3)))  # Outputs: Intruder detected! Math doesn't match reversed logic
