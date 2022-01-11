import numpy as np
import sys
import random
import math

from scipy.stats import ortho_group

n = 10
if len(sys.argv) == 2:
    n = int(sys.argv[1])
print("n=%d"%n)


A = np.random.rand(n, n) * 10
for i in range(n):
    for j in range(n):
        if i>j:
            A[i, j] = A[j, i]

Q, R = np.linalg.qr(A)

# print(Q.dot(R) - A)

Q_file = "input/Q.txt"
R_file = "input/R.txt"
A_file = "input/A.txt"

# np.savetxt(Q_file, Q, fmt="%f")
# np.savetxt(R_file, R, fmt="%f")
np.savetxt(A_file, A, fmt="%f")

with open(A_file, "r") as f:
    s = f.read();

with open(A_file, "w") as f:
    head = "%d %d\n"%(n, n)
    f.write(head + s)