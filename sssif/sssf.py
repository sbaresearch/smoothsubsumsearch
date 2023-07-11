"""Filtered Smooth Subsum Search (SSSf) algorithm as described in Section 3.3 
of 'Smooth Subsum Search: A Heuristic for Practical Integer Factorization'"""

import sys
import gmpy2
import time as tm
from collections import Counter
from SSSIF.mstep import matrix_factor
from math import prod, ceil, log, gcd
from random import seed, choices, randint
from sympy import sieve, isprime
from sympy.ntheory.residue_ntheory import legendre_symbol, sqrt_mod



def SSS(N, prop, digred):
    """Calling SSSf with two arguments prop and digred
    prop: Proportion of factor base used in filter
    digred: Determines bound for filtered candidates"""
    
    def pol(j):
        
        return (j + b)**2 - a * N


    def precomp():

        roots=[[] for _ in range(len(fbase))]
        for i in range(len(fbase)):
              p = fbase[i]
              roots[i].extend([(root-b) % p
                    for root in sqrt_mod(a*N, p, all_roots=True)])
              roots[i].sort()
        
        difflist=[coeffs[i]*(roots[i][1]-roots[i][0]) if i>0 else 0 
                 for i in range(len(plist))]
        
        return roots, difflist
    
    
    def products(xs):
        """Product tree computation as described in 
        https://facthacks.cr.yp.to/product.html"""
        result = [xs]
        while len(xs) > 1:
            xs = [prod(xs[i*2 : (i+1)*2]) for i in range((len(xs)+1)//2)]
            result.append(xs)
  
        return result


    def remainders(n, xlist):
        """Product tree computation as described in 
        https://facthacks.cr.yp.to/remainder.html"""
        x_tree = products(xlist)
        result = [n]
        for node in reversed(x_tree):
            result = [result[i//2] % node[i] for i in range(len(node))]
            
        return result
    
    
    def smooth_filter(xlist, val_list, zval):

        rem_list = remainders(zval, val_list)
        x_smooth, smooth_parts, filter_x, filter_v = [], [], [], []
        for i in range(len(rem_list)):
            y=rem_list[i]
            if y==0:
                x_smooth.append(xlist[i])
            else:
                n_smooth = val_list[i] // gcd(val_list[i], y)
                if n_smooth < filterbound:
                    if n_smooth <= last: 
                        smooth_parts.append([n_smooth, xlist[i]])
                    else: 
                        filter_x.append(xlist[i])
                        filter_v.append(n_smooth)
         
        return x_smooth, smooth_parts, filter_x, filter_v
    
    
    def smooth_batch(xlist, val_list, zval):

        rem_list = remainders(zval, val_list)
        x_smooth, smooth_parts = [], []
        for i in range(len(rem_list)):
            y=rem_list[i]
            if y==0:
                x_smooth.append(xlist[i])
            else:
                n_smooth = val_list[i] // gcd(val_list[i], y)
                if n_smooth < partialbound:
                    if n_smooth <= last: x_smooth.append(xlist[i])
                    else: smooth_parts.append([n_smooth, xlist[i]])
         
        return x_smooth, smooth_parts
    
       
    def search(output, partial, length):
        
        indset = set(choices(range(len(plist)), k=length))
        primes = [1] + [plist[i] for i in indset]
        indlist = [-1] + [i for i in indset]
        
        M = prod(primes)
        mvals = [M//p for p in primes] 
        xval = sum(coeffs[i]*roots[i][0] for i in indset) % M
        inv = [pow(M,-1,fbase[i]) if i>=len(plist) else 0 
               for i in range(len(fbase))]
        
        for i in indset:
            if i==0: continue
            xval = (xval+difflist[i]) % M  
            
            bvals, kvals = [0]*len(fbase), []
            for v in range(len(plist),len(fbase)):
                p, r = fbase[v], roots[v]
                b0 = (r[0] - xval)*inv[v] % p
                b1 = (r[1] - xval)*inv[v] % p
                bvals[v] = [b0, b1]
                kvals.extend([b0-p, b0, b1-p, b1])
             
            xlist, val_list = [], []
            for j in range(len(mvals)):
                if indlist[j]==i: continue
                m, pri = mvals[j], primes[j]
                
                if j>0:
                    kvals=[]
                    for v in range(len(plist), len(fbase)):
                        p, b = fbase[v], bvals[v]
                        r0 = b[0]*pri % p
                        r1 = b[1]*pri % p
                        kvals.extend([r0-p, r0, r1-p, r1])
                            
                for k,v in Counter(kvals).items():
                    if v>2:
                        arg=xval+k*m
                        xlist.append(arg)
                        val_list.append(abs(pol(arg))//m)
        
            sm, smp, fx, fv = smooth_filter(xlist, val_list, mval)
            res = smooth_batch(fx, fv, lval)
            sm.extend(res[0])
            smp.extend(res[1])
            
            output.update(sm)
            partial.extend(smp)
            
     
    """########################INITIALIZATION########################"""
    
        
    start = tm.time()
    
    if isprime(N): sys.exit("N is prime")
        
    if N < 2: sys.exit("N must be greater than 1")
    for r in range(2, N.bit_length()+1):
        x = round(log(N,r),10)
        if x == round(x): sys.exit("N is power of "+str(r))
        
    dig = len(str(N))  
    if dig <= 18: m = 120
    elif dig <= 25: m = 300
    elif dig <= 34: m = 400
    elif dig <= 36: m = 600
    elif dig <= 38: m = 800
    elif dig <= 40: m = 1000
    elif dig <= 42: m = 1200
    elif dig <= 44: m = 1400
    elif dig <= 48: m = 2000
    elif dig <= 52: m = 2400
    elif dig <= 56: m = 4000
    elif dig <= 60: m = 8000
    elif dig <= 66: m = 12000
    elif dig <= 70: m = 20000
    elif dig <= 74: m = 20000   
    elif dig <= 80: m = 60000
    elif dig <= 88: m = 100000
    elif dig <= 94: m = 120000
    else:           m = 200000
     
    a=1
    b=int(ceil(gmpy2.sqrt(a*N)))
    
    """Computing factor base"""
    sieve._reset()
    sieve.extend_to_no(m)
    fbase=list(sieve._list)[:m]
    plist=fbase[:m//5] 
    for p in fbase:
        if N % p == 0: sys.exit("Small prime divisor found: "+str(p))
    np_ind=[]
    for i in range(1,len(fbase)):
        if legendre_symbol(a*N,fbase[i])!=1: np_ind.append(i)    
    fbase=[fbase[i] for i in range(len(fbase)) if i not in np_ind]
    plist=[plist[i] for i in range(len(plist)) if i not in np_ind]   
    last = fbase[-1]
    
    m_tree=products(fbase[:len(fbase)//prop])
    l_tree=products(fbase[len(fbase)//prop:])
    n_tree=products(plist)
    mval = m_tree[-1][0]
    lval = l_tree[-1][0]
    nval = n_tree[-1][0]
    
    partialbound = 128*last
    filterbound = 10**(dig//2-digred)
    
    """Raising exponents in zvals for powersmooth test"""
    for p in fbase[:len(fbase)//prop]: 
        mval = mval * (p**(max(1,int(log(2**15,p)))-1))
    for p in fbase[len(fbase)//prop:]: 
        lval = lval * (p**(max(1,int(log(2**15,p)))-1))

    """Computing global coefficients for chinese remaindering"""      
    coeffs=[(pow(nval//plist[i],-1,plist[i])*nval//plist[i]) 
            for i in range(len(plist))]
    
    roots, difflist = precomp()
    print("N:", str(N)+" ("+str(len(str(N)))+" digits)")
    print("# Primes in FB: "+str(len(fbase)))
    print(" ")
    print("------------")
    print(" ")
    

    """########################MAIN LOOP########################"""
    
        
    print("*** Phase 1: Smooth Subsum Search ***")
        
    partial=[]
    lp_lst = []
    output = set()
    for run in range(10**100):
            
        search(output, partial, 7)
        pr=len(output)
        
        """Updating partial relations"""
        for pa in partial:
            l_p=[lst[0] for lst in lp_lst]
            if pa[0] in l_p:
                plst=lp_lst[l_p.index(pa[0])]
                if pa[1] not in plst: plst.append(pa[1])
            else:
                lp_lst.append(pa) 
        partial=[]
        """Keeping track of number of suitable partial relations"""
        lp=0
        for tupl in [lst for lst in lp_lst if len(lst)>2]:
            lp += (len(tupl)-1)//2
                 
        t=round(tm.time()-start,2)
        frac=(pr+lp+1)/(len(fbase)+10)
        print('\b'*512 +"Relations: "+str(pr+lp)+"/"+str(len(fbase)+10)+
              " ("+str("%0.2f" % min(100,round(frac*100,2)))+"%)"      
        +", running for "+str("%0.0f" % t)+"s   ", end='', flush=True)
         
        if pr + lp >= len(fbase) + 10: break
        
    large_primes=[lst for lst in lp_lst if len(lst)>2]
    s_nums=[pol(x) for x in output]
    xlist=[x + b for x in output]
     
    """Computing full relations from partial relations"""
    for tupl in large_primes:
        g=gcd(tupl[0],N)
        if g==1:
            i=1
            pairs=[]
            while i <= len(tupl[1:])-1:
                pairs.append([tupl[i],tupl[i+1]])
                i+=2  
            for pa in pairs:
                xlist.append((pa[0]+b)*(pa[1]+b) * pow(tupl[0],-1,N))
                s_nums.append(int(pol(pa[0]) * pol(pa[1]) // (tupl[0]**2)))
        else:
            sys.exit("Proper factors found: " + str(g) + " | " + str(N//g))
        
    print(" ")
    print("SSS finished: "+str(pr+lp)+" relations found," 
          +" "+str(lp)+ " of which from partial relations...")
    print("SSS time: "+str(round(tm.time()-start, 5)))
    print(" ")
    print("------------")
    print(" ")
        
    matrix_factor(N, xlist, s_nums, fbase)
        
    print("Total Time: "+str(round(tm.time()-start, 5)))



if __name__ == "__main__":
    
    seed(1)
    digits=75
    while True:
        
        if digits % 2 == 1:
            l=10**(digits//2-1)
            u=5*10**(digits//2)
        if digits % 2 == 0:
            l=10**(digits//2-1)
            u=5*10**(digits//2-1)
            
        p = randint(l,u)
        while not isprime(p): p = randint(l,u)
        q = randint(5*l,5*u)
        while not isprime(q): q = randint(5*l,5*u)    
        N=p*q
        dig = len(str(N))   
        if dig==digits: break
    
    SSS(N, 10, 5)



    