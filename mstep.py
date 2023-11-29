"""Linear algebra step for finding dependencies in the exponent vectors of the 
smooth relations. Implementation of Gaussian elimination is as in the primefac 
package: https://pypi.org/project/primefac/ See also: Çetin K. Koç and 
Sarath N. Arachchige. 'A Fast Algorithm for Gaussian Elimination over GF(2) 
and its Implementation on the GAPP.' Journal of Parallel and Distributed 
Computing 13.1 (1991): 118-122]"""


from time import time
from math import gcd
from gmpy2 import isqrt


def matrix_factor(N, xlist, smooth_nums, factor_base): 
        """xlist: List of the solutions x for which f(x) is smooth"""
        """smooth_nums: List of smooth values f(x) for x in xlist"""
    
        def factorize(n, factor_base):
        
            factors = []
            if n < 0:
                factors.append(-1)
            for p in factor_base:
                if p == -1: pass
                else:
                    while n % p == 0:
                        n //= p
                        factors.append(p)
            return factors
        
        factor_base.insert(0,-1)
        
        smooth_relations=[]
        for i in range(len(smooth_nums)):
            n=smooth_nums[i]
            n_factors = factorize(n,factor_base)
            divisors=[]
            for j in range(len(factor_base)):
                if factor_base[j] in n_factors:
                    divisors.append((j, n_factors.count(factor_base[j])))
            
            smooth_relations.append((xlist[i], n, divisors))
 
        relcount = len(smooth_relations)
        nf = len(factor_base)
    
        print("*** Phase 2: Linear Algebra ***")
        print("Building matrix...")
        M = [0] * nf
        mask = 1
        for sr in smooth_relations:
            for (j,exp) in sr[2]:
                if exp % 2: M[j] += mask
            mask <<= 1

        print("Gauss elimination...")
        gaussstart = time()
        row_is_marked = bytearray([False]) * relcount
        pivots = [-1] * nf
        for j in range(nf):
            M_j = M[j]
            i = (M_j & (-M_j)).bit_length() - 1 
            if i > -1:
                pivots[j] = i
                row_is_marked[i] = True
                for k in range(nf):
                    if (M[k] >> i) & 1 and k != j: 
                        M[k] ^= M_j

        print("Gaussian elimination time: %f seconds" % (time()-gaussstart))
        attempts = 0
        for i in range(relcount):
            if not row_is_marked[i]:
                square_indices = [i]
                for j in range(nf):
                    if (M[j] >> i) & 1:  
                        square_indices.append(pivots[j])
                attempts += 1
                sqrt1, sqrt2 = 1, 1
                for idx in square_indices:
                    sqrt1 *= smooth_relations[idx][0]
                    sqrt2 *= smooth_relations[idx][1]
                sqrt2 = isqrt(sqrt2)
                assert (sqrt1 * sqrt1) % N == (sqrt2 * sqrt2) % N
                factor = gcd(sqrt1 - sqrt2, N)
                if 1 != factor != N: 
                    print("Proper factors found: "
                          +str(factor) + " | " + str(N//factor))
                    return(factor)
        
        print("No proper factor found")