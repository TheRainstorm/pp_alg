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
B = np.random.rand(n, n) * 10
C = A.dot(B)

A_file = "input/A.txt"
B_file = "input/B.txt"
C_file = "input/C.txt"

np.savetxt(C_file, C, fmt="%f")

np.savetxt(A_file, A, fmt="%f")

with open(A_file, "r") as f:
    s = f.read();

with open(A_file, "w") as f:
    head = "%d %d\n"%(n, n)
    f.write(head + s)


np.savetxt(B_file, B, fmt="%f")

with open(B_file, "r") as f:
    s = f.read();

with open(B_file, "w") as f:
    head = "%d %d\n"%(n, n)
    f.write(head + s)