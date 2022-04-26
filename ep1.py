import numpy as np
import math
import sys
import time

def gerador_sistema_teste(n):
    i = np.multiply(np.ones((1,n)), np.resize(np.arange(1,n+1),(1,n)))
    
    a = (2*i -1)/(4*i)
    a[0][n-1] = (2*n - 1)/(2*n)
    
    b = np.ones((1,n)) * 2
    
    c = 1 - a

    d = np.cos((2*math.pi*(np.multiply(i,i)))/(n**2))

    return a,b,c,d

def LU_decomp_tridiag(a,b,c):
    n = len(b[0])
    u = np.zeros((1,n))
    l = np.zeros((1,n))
    # decompoe A em LU
    u[0][0] = b[0][0]
    l[0][0] = 0
    for i in range(1,n):
        l[0][i] = a[0][i]/u[0][i-1]
        u[0][i] = b[0][i]-l[0][i]*c[0][i-1]
    return u, l

def resolve_sistema_nao_ciclico(u,l,c,d):
    # resolve Ly=d
    n = len(u[0])
    y = np.zeros((1,n))
    x = np.zeros((1,n))
    y[0][0]=d[0][0]
    for i in range(1,n):
        y[0][i] = d[0][i]-l[0][i]*y[0][i-1]
    # resolve Ux = y
    x[0][n-1] = y[0][n-1]/u[0][n-1]
    for i in range(n-2,-1,-1):
        x[0][i]=(y[0][i]-c[0][i]*x[0][i+1])/u[0][i]
    return x
    
def resolve_sistema_ciclico(a,b,c,d):
    v = np.zeros((1,len(a[0])-1)) #is the size correct?
    v[0][0] = a[0][0]
    v[0][len(v[0])-1] = c[0][len(c[0])-2]
    v = v.T

    w = np.zeros((1,len(a[0])-1)) #is the size correct?
    w[0][0] = c[0][len(c[0])-1]
    w[0][len(w[0])-1] = a[0][len(a[0])-1]
    w = w.T

    d_til = d.T[0:len(d[0])-1]

    #Composicao da submatriz T
    a[0][0] = 0
    a = np.resize(a[0][0:len(a[0])-1], (1,len(a[0])-1))
    b = np.resize(b[0][0:len(b[0])-1], (1,len(b[0])-1))
    c = np.resize(c[0][0:len(c[0])-1], (1,len(c[0])-1))

    #Resolucao de T*y_til = d_til
    u, l = LU_decomp_tridiag(a,b,c)
    y_til = resolve_sistema_nao_ciclico(u,l,c,d_til.T) 
    
    #Resolucao de T*z_til = v
    z_til = resolve_sistema_nao_ciclico(u,l,c,v.T)
    print(z_til)

if __name__ == '__main__':
    #a = [0,-1,-1,-1]
    #b = [2, 2, 2, 2]
    #c = [-1, -1, -1, -1]
    #d = [1,0,0,1]

    a,b,c,d = gerador_sistema_teste(4)
    print("a:" + str(a))
    print("b:" + str(b))
    print("c:" + str(c))
    print("d:" + str(d))
    print("ahoy")
    
    start_time = time.time()
    u,l  = LU_decomp_tridiag(a,b,c)
    print(u)
    print(l)
    solucao = resolve_sistema_nao_ciclico(u,l,c,d)
    print(solucao)
    print(time.time() - start_time)
    resolve_sistema_ciclico(a,b,c,d)
