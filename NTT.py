#########################################################
#                                                       #
# File Name: NTT.py                                     #
#                                                       #
# Purpose: Implement number theoretic transform         #
#                                                       #
# Creation Date: 2017/12/29                             #
#                                                       #
# Last Modified: 2018/02/13                             #
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

    def bitReverse(self, num, len):
        """
        integer bit reverse
        input: num, bit length
        output: rev_num 
        example: input 6(110) output 3(011)
        complexity: O(len)
        """
        rev_num = 0
        for i in range(0, len):
            if (num >> i) & 1:
                rev_num |= 1 << (len - 1 - i)
        return rev_num

    def orderReverse(self, poly, N_bit):
        """docstring for order"""
        for i, coeff in enumerate(poly):
            rev_i = self.bitReverse(i, N_bit)
            if rev_i > i:
                coeff ^= poly[rev_i]
                poly[rev_i] ^= coeff
                coeff ^= poly[rev_i]
                poly[i] = coeff
        return poly

    # Basic method of NTT
    # The complexity is O(N log N)
    def ntt(self, poly, M, N, w):
        """number theoretic transform algorithm"""
        N_bit = N.bit_length() - 1
        rev_poly = self.orderReverse(poly, N_bit)
        for i in range(0, N_bit):
            points1, points2 = [], []
            for j in range(0, int(N / 2)):
                shift_bits = N_bit - 1 - i
                P = (j >> shift_bits) << shift_bits
                w_P = self.modExponent(w, P, M)
                odd = poly[2 * j + 1] * w_P
                even = poly[2 * j]
                points1.append((even + odd) % M)
                points2.append((even - odd) % M)
                # TODO: use barrett modular reduction
                points = points1 + points2
            if i != N_bit:
                poly = points
        return points

    # Basic metho of inverse NTT
    # The complexity is O(N log N)
    def intt(self, points, M, N, w):
        """inverse number theoretic transform algorithm"""
        inv_w = self.modInv(w, M)
        inv_N = self.modInv(N, M)
        p = []
        poly = self.ntt(points, M, N, inv_w)
        for i in range(0, N):
            poly[i] = (poly[i] * inv_N) % M
            # TODO: use barrett modular reduction
        return poly
