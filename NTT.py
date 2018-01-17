#########################################################
#                                                       #
# File Name: NTT.py                                     #
#                                                       #
# Purpose: generate primitive Nth root of unity         #
#                                                       #
# Creation Date: 2017/12/29                             #
#                                                       #
# Last Modified: 2017/01/10                             #
#                                                       #
# Create By: Jyun-Neng Ji                               #
#                                                       #
#########################################################
import math
import random
from sympy import isprime


class NTT:

    def isInteger(self, M):
        return type(M).__name__ == 'int'

    def isPrime(self, M):
        assert self.isInteger(M), 'Not an integer.'
        return isprime(M)

    # factorize algorithm
    # complexity is O(N^(1/2))
    # not used
    def factorize(self, N):
        factors = []
        for factor in range(1, int(math.sqrt(N) + 1)):
            if N % factor == 0:
                if factor != 1:
                    factors.append(factor)
                    factors.append(int(N / factor))
        factors.sort()
        return factors

    # modular expoential algorithm
    # complexity is O(log N)
    def modExponent(self, base, power, M):
        result = 1
        power = int(power)
        base = base % M
        while power > 0:
            if power & 1:
                result = (result * base) % M
            base = (base * base) % M
            power = power >> 1
        return result

    # calculate x^(-1) mod M
    def modInv(self, x, M):
        t, new_t, r, new_r = 0, 1, M, x

        while new_r != 0:
            quotient = int(r / new_r)
            tmp_new_t = new_t
            tmp_t = t
            tmp_new_r = new_r
            tmp_r = r
            t, new_t = new_t, (t - quotient * new_t)
            r, new_r = new_r, (r % new_r)
        if r > 1:
            return "x is not invertible."
        if t < 0:
            t = t + M
        return t

    # check if r^k = 1 (mod M), k<N
    def existSmallN(self, r, M, N):
        for k in range(2, N):
            if self.modExponent(r, k, M) == 1:
                return True
        return False

    # generate primitive nth root of unity
    # let a is primitive root i.e a^phi(M) = 1 (mod M) then
    # B = a^(phi(M)/N) (mod M) is Nth root of unity iff N|phi(M)
    def NthRootOfUnity(self, M, N):
        assert self.isPrime(M), 'Not a prime.'  # modulus should be a prime
        assert (M - 1) % N == 0, 'N cannot divide phi(M)'
        phi_M = M - 1
        while True:
            alpha = random.randrange(1, M)
            beta = self.modExponent(alpha, phi_M / N, M)
            # check if beta can be k th root of unity, k<N
            if not self.existSmallN(beta, M, N):
                return int(beta)

    # verify B^N = 1 (mod M)
    def isNthRootOfUnity(self, M, N, beta):
        return self.modExponent(beta, N, M) == 1

    # Basic method of NTT 
    # The complexity is O(N^2)
    def ntt(self, poly, M, N, w):
        """number theoretic transform algorithm"""
        points = []
        for k in range(0, N):
            point = 0
            for i in range(0, N):
                point += poly[i] * self.modExponent(w, i * k, M)
            point = point % M
            points.append(point)
        return points
    # TODO: reduce it's complexity

    # Basic metho of inverse NTT
    # The complexity is O(N^2)
    def intt(self, points, M, N, w):
        """inverse number theoretic transform algorithm"""
        poly = []
        inv_w = self.modInv(w, M)
        inv_N = self.modInv(N, M)
        for i in range(0, N):
            coeff = 0
            for k in range(0, N):
                coeff += points[k] * self.modExponent(inv_w, i * k, M)
            coeff *= inv_N
            coeff %= M
            poly.append(coeff)
        return poly
    # TODO: reduce it's complexity
