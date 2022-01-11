import numpy as np
import sys

n = 10
if len(sys.argv) == 2:
    n = int(sys.argv[1])
print("n=%d"%n)

L = np.random.rand(n, n)
U = np.random.rand(n, n)

for i in range(n):
    for j in range(n):
        if i > j:
            L[i, j] = 0
        if i < j:
            U[i, j] = 0

L_file = "input/L.txt"
U_file = "input/U.txt"
A_file = "input/A.txt"

np.savetxt(L_file, L, fmt="%f")
np.savetxt(U_file, U, fmt="%f")

A = np.matmul(L, U)
np.savetxt(A_file, A, fmt="%f")

with open(A_file, "r") as f:
    s = f.read();

with open(A_file, "w") as f:
    head = "%d %d\n"%(n, n)
    f.write(head + s)