import qutip as qt
import numpy as np

def quantum_kernel(x1, x2):
    state1 = qt.Qobj(np.cos(x1) * np.array([1, 0]) + np.sin(x1) * np.array([0, 1]))
    state2 = qt.Qobj(np.cos(x2) * np.array([1, 0]) + np.sin(x2) * np.array([0, 1]))
    return abs(state1.overlap(state2))**2

x = np.array([0.1, 0.2])
K = np.zeros((2,2))
for i in range(2):
    for j in range(2):
        K[i,j] = quantum_kernel(x[i], x[j])

print("Quantum Kernel Matrix:")
print(K)
