// numberSpeaker.js

const SMALL = [
  "zero","one","two","three","four","five","six","seven","eight","nine",
  "ten","eleven","twelve","thirteen","fourteen","fifteen","sixteen",
  "seventeen","eighteen","nineteen"
];

const TENS = [
  "", "", "twenty","thirty","forty","fifty","sixty","seventy","eighty","ninety"
];

// Extend this list as far as you want — JS can handle infinite scales
const SCALES = [
  "", "thousand", "million", "billion", "trillion",
  "quadrillion", "quintillion", "sextillion", "septillion",
  "octillion", "nonillion", "decillion", "undecillion",
  "duodecillion", "tredecillion", "quattuordecillion",
  "quindecillion", "sexdecillion", "septendecillion",
  "octodecillion", "novemdecillion", "vigintillion"
];

function threeDigitToWords(n) {
  n = Number(n);
  if (n === 0) return "";

  let words = [];

  const hundreds = Math.floor(n / 100);
  const remainder = n % 100;

  if (hundreds > 0) {
    words.push(SMALL[hundreds]);
    words.push("hundred");
  }

  if (remainder > 0) {
    if (remainder < 20) {
      words.push(SMALL[remainder]);
    } else {
      const tens = Math.floor(remainder / 10);
      const ones = remainder % 10;

      words.push(TENS[tens]);
      if (ones > 0) words.push(SMALL[ones]);
    }
  }

  return words.join(" ");
}

function numberToWords(num) {
  // Accept numbers or strings
  let str = String(num).trim();

  // Handle sign
  let sign = "";
  if (str.startsWith("-")) {
    sign = "negative ";
    str = str.slice(1);
  }

  // Remove leading zeros
  str = str.replace(/^0+/, "");
  if (str.length === 0) return "zero";

  // Split into 3-digit groups from the right
  let groups = [];
  while (str.length > 0) {
    groups.push(str.slice(-3));
    str = str.slice(0, -3);
  }

  let parts = [];

  for (let i = 0; i < groups.length; i++) {
    const groupVal = Number(groups[i]);
    if (groupVal === 0) continue;

    const chunk = threeDigitToWords(groupVal);
    const scale = SCALES[i] || `10^${i * 3}`;

    parts.push(`${chunk} ${scale}`.trim());
  }

  return (sign + parts.reverse().join(" ")).trim();
}

// Example usage:
console.log(numberToWords("9876543210123456789"));
console.log(numberToWords(42));
console.log(numberToWords("1000000000000000000000001"));
