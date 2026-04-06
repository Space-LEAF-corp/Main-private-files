# number_speaker.py

SMALL_NUMS = [
    "zero", "one", "two", "three", "four", "five", "six",
    "seven", "eight", "nine", "ten", "eleven", "twelve",
    "thirteen", "fourteen", "fifteen", "sixteen",
    "seventeen", "eighteen", "nineteen"
]

TENS = [
    "", "", "twenty", "thirty", "forty",
    "fifty", "sixty", "seventy", "eighty", "ninety"
]

# You can extend this list as far as you like
SCALES = [
    "", "thousand", "million", "billion", "trillion",
    "quadrillion", "quintillion", "sextillion", "septillion",
    "octillion", "nonillion", "decillion", "undecillion",
    "duodecillion", "tredecillion", "quattuordecillion",
    "quindecillion", "sexdecillion", "septendecillion",
    "octodecillion", "novemdecillion", "vigintillion"
]


def _three_digit_to_words(n: int) -> str:
    """Convert a number from 0–999 into words."""
    assert 0 <= n <= 999

    if n == 0:
        return ""

    words = []

    hundreds = n // 100
    remainder = n % 100

    if hundreds > 0:
        words.append(SMALL_NUMS[hundreds])
        words.append("hundred")

    if remainder > 0:
        if remainder < 20:
            words.append(SMALL_NUMS[remainder])
        else:
            tens = remainder // 10
            ones = remainder % 10
            words.append(TENS[tens])
            if ones > 0:
                words.append(SMALL_NUMS[ones])

    return " ".join(words)


def number_to_words(num) -> str:
    """
    Convert an arbitrarily large integer (int or string) to English words.
    Example: "9876543210123456789" ->
    "nine quintillion eight hundred seventy six quadrillion
     five hundred forty three trillion two hundred ten billion
     one hundred twenty three million four hundred fifty six thousand
     seven hundred eighty nine"
    """
    # Accept int or string
    if isinstance(num, int):
        num_str = str(num)
    else:
        num_str = str(num).strip()

    # Handle sign
    sign = ""
    if num_str.startswith("-"):
        sign = "negative "
        num_str = num_str[1:].lstrip()

    # Strip leading zeros
    num_str = num_str.lstrip("0")
    if num_str == "":
        return "zero"

    # Split into groups of 3 digits from the right
    groups = []
    while num_str:
        groups.append(num_str[-3:])
        num_str = num_str[:-3]
    # groups[0] is lowest (units), groups[1] is thousands, etc.

    words_chunks = []
    for idx, group_str in enumerate(groups):
        group_val = int(group_str)
        if group_val == 0:
            continue

        chunk_words = _three_digit_to_words(group_val)
        scale_word = SCALES[idx] if idx < len(SCALES) else f"10^{idx*3}"

        if scale_word:
            words_chunks.append(f"{chunk_words} {scale_word}".strip())
        else:
            words_chunks.append(chunk_words)

    # Highest scale first
    full = " ".join(reversed(words_chunks))
    return (sign + full).strip()


if __name__ == "__main__":
    tests = [
        "0",
        "7",
        "42",
        "100",
        "999",
        "1000",
        "123456789",
        "9876543210123456789",
        "-123456789012345678901234567890"
    ]

    for t in tests:
        print(t, "->", number_to_words(t))
