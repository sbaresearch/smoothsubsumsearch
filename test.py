import os
import time as tm

import numpy as np
from sssif import sss, sssf
from ssiqs import qs
from psiqs import siqs
from random import randint

from sympy import isprime
from sympy.ntheory.generate import prime

file_path = os.path.dirname(__file__)

runs=2
digits=70
all_algs=True

sss_t=[]
sssf_t=[]
siqs_t=[]
sysiqs_t=[]

n_list=[]
for i in range(runs):
    
    print(" ")
    print("***********TRIAL "+str(i+1)+"***********")
    print(" ")
    
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
    
    n_list.append(N)
    
    if digits < 60: 
        print("1: SSS")
        start=tm.time()
        sss.SSS(N)
        sss_t.append(tm.time()-start)
        
        print("\n\n\n")
    
    if digits > 55:
        print("1f: SSSF")
        start=tm.time()
        if 55<=digits<=60:
            sssf.SSS(N, 20, 8)
        if 61<=digits<=65:
            sssf.SSS(N, 22, 9)
        if 66<=digits<=70:
            sssf.SSS(N, 25, 10)
        if digits>70:
            sssf.SSS(N, 30, 12)
        sssf_t.append(tm.time()-start)
        
        print("\n\n\n")

    if all_algs:
        
        if digits < 60: 
            print("2: pSIQS")
            start=tm.time()
            p=siqs(N,verbose=True)
            print(p)
            siqs_t.append(tm.time()-start)
            
            print("\n\n\n")
        
        if   dig <= 34: nf, m = 200, 65536
        elif dig <= 36: nf, m = 300, 65536
        elif dig <= 38: nf, m = 400, 65536
        elif dig <= 40: nf, m = 500, 65536
        elif dig <= 42: nf, m = 600, 65536
        elif dig <= 44: nf, m = 700, 65536
        elif dig <= 48: nf, m = 1000, 65536
        elif dig <= 52: nf, m = 1200, 65536
        elif dig <= 56: nf, m = 2000, 65536 * 3
        elif dig <= 60: nf, m = 4000, 65536 * 3
        elif dig <= 66: nf, m = 6000, 65536 * 3
        elif dig <= 70: nf, m = 10000, 65536 * 3
        else:           nf, m = 30000, 65536 * 6
        
        print("3: sSIQS")
        start=tm.time()
        qs(N, prime(2*nf), m)
        sysiqs_t.append(tm.time()-start)
    
print(" ")
print("********************************Results********************************")    
if digits < 60:     
    print("SSS Mean: ", np.mean(sss_t))
    print("SSS STD: ", np.std(sss_t))
    print(" ")

if digits > 55:
    print("SSSF Mean: ", np.mean(sssf_t))
    print("SSSF STD: ", np.std(sssf_t))
    print(" ")


if all_algs:
    if digits < 60: 
        print("PSIQS Mean: ", np.mean(siqs_t))
        print("PSIQS STD: ", np.std(siqs_t))
        print(" ")
    print("SYSIQS Mean: ", np.mean(sysiqs_t))
    print("SYSIQS STD: ", np.std(sysiqs_t))
    
    table = os.path.join(file_path+'/results',
                    'timings_'+str(runs)+'_'+str(digits)+'.txt')
    
    with open(table, 'w') as w:
        if digits < 60: 
            w.write(
                "SSS Mean: "+str(np.mean(sss_t))+'\n'
                "SSS STD: "+str(np.std(sss_t))+'\n'
                '\n')
        if digits > 55:
            w.write(
            "SSSF Mean: "+str(np.mean(sssf_t))+'\n'
            "SSSF STD: "+str(np.std(sssf_t))+'\n'
            '\n')
        if digits < 60: 
            w.write("PSIQS Mean: "+str(np.mean(siqs_t))+'\n'
                    "PSIQS STD: "+str(np.std(siqs_t))+'\n'
                    '\n')
        w.write("SYSIQS Mean: "+str(np.mean(sysiqs_t))+'\n'
        "SYSIQS STD: "+str(np.std(sysiqs_t))+'\n')
        w.write("\n")
        w.write("***********************\n")
        w.write("List of N:\n")
        for n in n_list:
            w.write(str(n)+"\n")
        w.write("***********************\n")
        w.write("\n")
        if digits < 60: 
            w.write("SSS Timings: "+str(sss_t)+"\n")
        if digits > 55:
            w.write("SSSF Timings: "+str(sssf_t)+"\n")
        if digits < 60: 
            w.write("pSIQS Timings: "+str(siqs_t)+"\n")
        w.write("sSIQS Timings: "+str(sysiqs_t)+"\n")
        w.write("\n")
        



























