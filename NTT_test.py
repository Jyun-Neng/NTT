from NTT import NTT

ntt = NTT()
M = 2013265921
N = 2048
print("Modulus is %d" %M)
print("The number of FFT points is %d" %N)
r = ntt.NthRootOfUnity(M, N)
print ("The %dth root of unity is %d" % (N, r))
if ntt.isNthRootOfUnity(M, N, r):
    print("%d to power of %d is congruent to 1 modulo %d" % (r, N, M))
else:
    print("%d to power of %d is not congruent to 1 modulo %d" % (r, N, M))
