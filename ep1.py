import numpy as np
import math
import sys
import time

def gerador_sistema_teste(n):
    #criação de um vetor crescente [1,2,3,..., n]
    i = np.resize(np.arange(1,n+1),(1,n))
    
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
    
def resolve_sistema_ciclico(a,b,c,d,n):
    v = np.zeros((1,len(a[0])-1)) #is the size correct?
    v[0][0] = a[0][0]
    v[0][len(v[0])-1] = c[0][len(c[0])-2]
    v = v.T

    w = np.zeros((1,len(a[0])-1)) #is the size correct?
    w[0][0] = c[0][len(c[0])-1]
    w[0][len(w[0])-1] = a[0][len(a[0])-1]
    w = w.T

    d_til = d.T[0:len(d[0])-1]
    print(d_til)

    #Composicao da submatriz T
    a_T = np.resize(a[0][0:len(a[0])-1], (1,len(a[0])-1))
    a_T[0][0] = 0
    b_T = np.resize(b[0][0:len(b[0])-1], (1,len(b[0])-1))
    c_T = np.resize(c[0][0:len(c[0])-1], (1,len(c[0])-1))
    c_T[0][len(c_T[0])-1] = 0

    #Resolucao de T*y_til = d_til
    u, l = LU_decomp_tridiag(a_T,b_T,c_T)
    y_til = resolve_sistema_nao_ciclico(u,l,c_T,d_til.T) 
    

    #Resolucao de T*z_til = v
    z_til = resolve_sistema_nao_ciclico(u,l,c_T,v.T)
    
    #Gera a  matriz T inteira para teste e calcula o residuo dos sistemas nao ciclicos
    T = gera_matriz_tridiagonal_n_cicl(len(a[0])-1)
    print("residuo de Ty = d: ")
    print(d_til-np.dot(T,y_til.T))
    print("residuo de Tz = v")
    print(v-np.dot(T,((z_til).T)))
    
    #Determinacao de xn
    xn = d[0][n-1] - (c[0][n-1])*(y_til[0][0]) - (a[0][n-1])*(y_til[0][n-2])
    xn = xn/(b[0][n-1] - (c[0][n-1])*(z_til[0][0]) - (a[0][n-1])*(z_til[0][n-2]))
    #Determinacao de x_til
    x_til = y_til - xn*z_til
    return np.resize(np.insert(x_til,len(x_til[0]),xn),(1,len(x_til[0])+1)).T

def gera_matriz_tridiagonal(n):
    a = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i==j+1 and i !=n-1:
                a[i][j] = (2*(i+1)-1)/(4*(i+1))
            if i == j:
                a[i][j] = 2
            if j == i+1:
                a[i][j] = 1-(2*(i+1)-1)/(4*(i+1))
    a[n-1][n-2] = (2*n-1)/(2*n)
    a[0][n-1] = 1/4
    a[n-1][0] = 1-(2*n-1)/(2*n)
    return a

def gera_matriz_tridiagonal_n_cicl(n):
    a = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i==j+1:
                a[i][j] = (2*(i+1)-1)/(4*(i+1))
            if i == j:
                a[i][j] = 2
            if j == i+1:
                a[i][j] = 1-(2*(i+1)-1)/(4*(i+1))
    return a

if __name__ == '__main__':
    #a = [0,-1,-1,-1]
    #b = [2, 2, 2, 2]
    #c = [-1, -1, -1, -1]
    #d = [1,0,0,1]
    n = 20
    a,b,c,d = gerador_sistema_teste(20)
    print("a:" + str(a))
    print("b:" + str(b))
    print("c:" + str(c))
    print("d:" + str(d))
    
    start_time = time.time()
    u,l  = LU_decomp_tridiag(a,b,c)
    x = resolve_sistema_ciclico(a,b,c,d,n)
    print("Resultado: ")
    print(x)
    print("Tempo decorrido: ", time.time() - start_time)
    A = gera_matriz_tridiagonal(20)
    r = np.dot(A,x) - d.T
    print("residuo final: ")
    print(r)
    x_barra  = np.array([[0.33189941,  0.33337847,  0.33089203,  0.32456897,  0.31054214,  0.28498039, 0.24375753,  0.1834913, 0.10274418,  0.00360625, -0.10669712, -0.21472831, -0.30113601, -0.34331333, -0.32095628, -0.22457847, -0.06361961,  0.12491935,0.29035857,  0.34417587]])
    e = x - (x_barra).T
    print((x_barra).T)
    print(e)
    print(r+A@e)
    


    #### TO DO ####
    # Documentação
        # Comentar cada coisinha do codigo
        # README
        # Relatorio
    # Ultimo resultado está errado? -> resíduo anormalmente grande
