import sys
import random

if len(sys.argv) < 2:
    print("this function takes an ergumet")
    exit()

repeats = 1
if len(sys.argv) > 2:
    repeats = int(sys.argv[2])

N = int(sys.argv[1])
A = [(a,b) for a in [-0.5, 0.0, 0.5, 1.0] for b in [-1.0,-0.5,0.0,0.5,1.0]]

print(N*repeats)
l = []
for i in range(N):
    a = A[random.randint(0,19)]
    l.append(a[0])
    l.append(a[1])
for j in range(repeats):
    for i in range(2*N):
        print(l[i])
