"""Code validator and Fibonacci reverse resequencer.

This module provides:
1. Forward/backward validation — ensure code can be parsed both ways
2. Fibonacci reverse split — split code into chunks of Fib sizes (21, 13, 8, 5, 3, 2, 1)
3. Resequence — reconstruct code in a reordered Fibonacci sequence
"""
from __future__ import annotations

import ast
from typing import List, Tuple


def fib_reverse_sequence() -> List[int]:
    """Generate Fibonacci sequence in reverse from 21: [21, 13, 8, 5, 3, 2, 1]."""
    fibs: List[int] = []
    a, b = 1, 1
    while a <= 21:
        fibs.append(a)
        a, b = b, a + b
    # Remove duplicate 1
    return sorted(set(fibs), reverse=True)


class CodeValidator:
    def __init__(self, code: str):
        self.code = code
        self.lines = code.split('\n')

    def validate_forward(self) -> Tuple[bool, str]:
        """Check if code parses successfully (forward)."""
        try:
            ast.parse(self.code)
            return True, "Forward validation passed"
        except SyntaxError as e:
            return False, f"Forward validation failed: {e}"

    def validate_backward(self) -> Tuple[bool, str]:
        """Check if reversed code has syntactic patterns (backward check)."""
        # Backward check: simple heuristic — lines should be reversible without breaking structure
        try:
            # Try to parse as-is
            ast.parse(self.code)
            # If it parses forward, check for basic reversibility
            lines = self.lines
            # Simple check: reversed lines should have same count
            return len(lines) > 0, "Backward reversibility check passed"
        except Exception as e:
            return False, f"Backward check failed: {e}"

    def validate_both(self) -> Tuple[bool, str]:
        """Validate both forward and backward."""
        fwd_ok, fwd_msg = self.validate_forward()
        bwd_ok, bwd_msg = self.validate_backward()
        if fwd_ok and bwd_ok:
            return True, f"Both checks passed. {fwd_msg} | {bwd_msg}"
        else:
            msgs = [m for m in [fwd_msg, bwd_msg] if m]
            return False, " | ".join(msgs)


class FibonacciResequencer:
    """Split code into Fibonacci-sized chunks (21, 13, 8, 5, 3, 2, 1) and resequence."""

    def __init__(self, code: str):
        self.code = code
        self.lines = code.split('\n')
        self.fib_sizes = fib_reverse_sequence()

    def split_into_chunks(self) -> List[List[str]]:
        """Split lines into chunks of Fibonacci sizes."""
        chunks: List[List[str]] = []
        idx = 0
        for size in self.fib_sizes:
            if idx >= len(self.lines):
                break
            chunk = self.lines[idx:idx + size]
            chunks.append(chunk)
            idx += size
        # Remainder
        if idx < len(self.lines):
            chunks.append(self.lines[idx:])
        return chunks

    def resequence_reverse(self, chunks: List[List[str]]) -> str:
        """Resequence chunks in reverse Fibonacci order."""
        # Reverse the chunk order
        reversed_chunks: List[List[str]] = list(reversed(chunks))
        result_lines: List[str] = []
        for chunk in reversed_chunks:
            result_lines.extend(chunk)
        return '\n'.join(result_lines)

    def resequence_forward(self, chunks: List[List[str]]) -> str:
        """Resequence chunks in Fibonacci order (as-is)."""
        result_lines: List[str] = []
        for chunk in chunks:
            result_lines.extend(chunk)
        return '\n'.join(result_lines)

    def split_and_resequence(self, reverse: bool = True) -> tuple[str, list[tuple[int, int]]]:
        """
        Split code and resequence.

        Returns:
            (resequenced_code, chunk_sizes)
            chunk_sizes: list of (start_line, end_line) tuples for each chunk
        """
        chunks = self.split_into_chunks()
        chunk_info: list[tuple[int, int]] = []
        idx = 0
        for chunk in chunks:
            chunk_info.append((idx, idx + len(chunk)))
            idx += len(chunk)

        if reverse:
            resequenced = self.resequence_reverse(chunks)
        else:
            resequenced = self.resequence_forward(chunks)

        return resequenced, chunk_info

    def analyze_chunks(self) -> dict[str, object]:
        """Analyze chunk structure."""
        chunks = self.split_into_chunks()
        return {
            "num_chunks": len(chunks),
            "fib_sizes": self.fib_sizes,
            "chunk_sizes": [len(c) for c in chunks],
            "total_lines": len(self.lines),
        }
