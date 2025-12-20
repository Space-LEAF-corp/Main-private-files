# day1_cipher_matrix.py
# Cryptographic matrix cipher (Day-1)
# Master strand: "CTTTTCTATCGA" duplicated to 24 bases

def reshape_to_matrix(seq: str, rows: int, cols: int) -> list[list[str]]:
    return [list(seq[i*cols:(i+1)*cols]) for i in range(rows)]

def flatten_matrix(matrix: list[list[str]]) -> str:
    return ''.join(''.join(row) for row in matrix)

def rotate_column_down(matrix: list[list[str]], col: int, steps: int = 1) -> None:
    col_vals: list[str] = [row[col] for row in matrix]
    steps = steps % len(col_vals)
    rotated: list[str] = col_vals[-steps:] + col_vals[:-steps]
    for i in range(len(matrix)):
        matrix[i][col] = rotated[i]

def rotate_column_up(matrix: list[list[str]], col: int, steps: int = 1) -> None:
    rotate_column_down(matrix, col, -steps)

from typing import Iterable

def forward_cipher(seq: str, mask_cols: Iterable[int], k: int, g1: int, g2: int) -> str:
    rows, cols = 4, 6
    matrix = reshape_to_matrix(seq, rows, cols)
    # Column stomp
    for c in mask_cols:
        rotate_column_down(matrix, c, k)
    # Global stomps
    for c in range(cols):
        rotate_column_down(matrix, c, g1)
        rotate_column_down(matrix, c, g2)
    return flatten_matrix(matrix)

def inverse_cipher(cipher_seq: str, mask_cols: Iterable[int], k: int, g1: int, g2: int) -> str:
    rows, cols = 4, 6
    matrix = reshape_to_matrix(cipher_seq, rows, cols)
    # Global lifts
    for c in range(cols):
        rotate_column_up(matrix, c, g2)
        rotate_column_up(matrix, c, g1)
    # Column lifts
    for c in mask_cols:
        rotate_column_up(matrix, c, k)
    return flatten_matrix(matrix)

if __name__ == "__main__":
    master = "CTTTTCTATCGA" * 2  # 24 bases
    mask_cols = [0, 2, 4]        # stomp columns (0-indexed)
    k, g1, g2 = 1, 1, 1          # stomp parameters

    cipher = forward_cipher(master, mask_cols, k, g1, g2)
    print("Day-1 Cipher View:", cipher)

    restored = inverse_cipher(cipher, mask_cols, k, g1, g2)
    print("Restored Master:", restored)
