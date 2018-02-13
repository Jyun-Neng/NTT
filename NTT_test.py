import numpy as np
from NTT import NTT

ntt = NTT()
M = 2013265921
N = 2048
r = ntt.NthRootOfUnity(M, N)
print("Modulus is %d" % M)
print("The number of FFT points is %d" % N)
print("The %dth root of unity is %d" % (N, r))
if ntt.isNthRootOfUnity(M, N, r):
    print("%d to power of %d is congruent to 1 modulo %d" % (r, N, M))
else:
    print("%d to power of %d is not congruent to 1 modulo %d" % (r, N, M))


M = 2013265921
poly = [1, 2, 3, 4]     # 4x^3+3x^2+2x+1
print("Modulus : %d" % M)
print("Polynomial : ", poly)
N = len(poly)
w = ntt.NthRootOfUnity(M, N)
ntt_poly = ntt.ntt(poly, M, N, w)
intt_poly = ntt.intt(ntt_poly, M, N, w)
print("Polynomial degree : %d" % (N - 1))
print("Primitive %dth root of unity : %d" % (N, w))
print("NTT(poly) = ", ntt_poly)
print("Inverse NTT : ", intt_poly)
