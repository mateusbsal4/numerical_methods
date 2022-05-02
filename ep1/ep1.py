import numpy as np
import math
import sys
import time

#Mateus Bonelli Salomão 11914789
#Pedro Pimentel Fuoco 11804313

def gerador_sistema_teste(n): #Essa funcao gera o sistema teste proposto no enunciado do exercicio programa, retornando a,b,c e d de uma matrix n por n, sendo n a entrada.
    
    #criaçao de um vetor crescente [1,2,3,..., n]
    i = np.resize(np.arange(1,n+1),(1,n))
    
    a = (2*i -1)/(4*i)
    a[0][n-1] = (2*n - 1)/(2*n)
    
    b = np.ones((1,n)) * 2
    
    c = 1 - a

    d = np.cos((2*math.pi*(np.multiply(i,i)))/(n**2))

    return a,b,c,d

def LU_decomp_tridiag(a,b,c): #Essa funcao recebe uma matriz tridiagonal em sua entrada, nas variaveis a,b e c, e retorna u e l, coeficientes das matrizes de decomposicao LU
    n = len(b[0])
    u = np.zeros((1,n))
    l = np.zeros((1,n))
    # Seguimos o algoritmo proposto no enunciado do EP1
    u[0][0] = b[0][0]
    l[0][0] = 0
    for i in range(1,n):
        l[0][i] = a[0][i]/u[0][i-1]
        u[0][i] = b[0][i]-l[0][i]*c[0][i-1]
    return u, l

def resolve_sistema_nao_ciclico(u,l,c,d): #Resolve para uma matriz  tridiagonal nao ciclica decomposta em LU. as entradas são u,l,c e d. devolve o vetor x, resolucao do sistema linear
    n = len(u[0])
    y = np.zeros((1,n))
    x = np.zeros((1,n))
    y[0][0]=d[0][0]
    for i in range(1,n):
        y[0][i] = d[0][i]-l[0][i]*y[0][i-1]
    x[0][n-1] = y[0][n-1]/u[0][n-1]
    for i in range(n-2,-1,-1):
        x[0][i]=(y[0][i]-c[0][i]*x[0][i+1])/u[0][i]
    return x
    
def resolve_sistema_ciclico(a,b,c,d,n): #Resolve Ax=d para matriz A tridiagonal ciclica. as entradas sao a,b,c,d e n (tamanho da matriz). devolve o vetor x, resolucao do sistema linear
    v = np.zeros((1,len(a[0])-1))
    v[0][0] = a[0][0]
    v[0][len(v[0])-1] = c[0][len(c[0])-2]
    v = v.T


    d_til = d.T[0:len(d[0])-1]

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
    
    #print("residuo de Ty_til = d_til: ")
    #print(d_til-np.dot(T,y_til.T))
    #print("residuo de Tz_til = v")
    #print(v-np.dot(T,((z_til).T)))
    
    #Gera a  matriz T inteira para teste e calcula o residuo dos sistemas nao ciclicos
    T = gera_matriz_tridiagonal_n_cicl(len(a[0])-1)
    
    
    #Determinacao de xn
    xn = d[0][n-1] - (c[0][n-1])*(y_til[0][0]) - (a[0][n-1])*(y_til[0][n-2])
    xn = xn/(b[0][n-1] - (c[0][n-1])*(z_til[0][0]) - (a[0][n-1])*(z_til[0][n-2]))
    
    #Determinacao de x_til
    x_til = y_til - xn*z_til
    
    #Retorna x_til com xn adicionado no final do vetor
    return np.resize(np.insert(x_til,len(x_til[0]),xn),(1,len(x_til[0])+1)).T

def gera_matriz_tridiagonal(n): #Retorna matriz tridiagonal nao ciclica de teste "a" de tamanho n por n. Essa matriz servira para as funcoes de debug do codigo
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

def gera_matriz_tridiagonal_n_cicl(n): #Retorna matriz tridiagonal ciclica de teste "a" de tamanho n por n. Essa matriz servira para as funcoes de debug do codigo
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
    modo = int(input("Digite o número do teste que você quer rodar, sendo: \n Teste 1: achar decomposiçao LU da matriz tridiagonal T \n Teste 2: Resolver sistemas nao ciclicos Ty_til = d_til e Tz_til=v  \n Teste 3: resolver sistema ciclico Ax = d \n"))
    if modo != 1 and modo !=2 and modo!= 3:
        print("Entrada invalida")
        sys.exit()
    n = int(input("Digite n: "))
    assert n > 0 and type(n) == int
    
   
    
    a,b,c,d = gerador_sistema_teste(n)
    
    
    a_T = np.resize(a[0][0:len(a[0])-1], (1,len(a[0])-1))
    a_T[0][0] = 0
    b_T = np.resize(b[0][0:len(b[0])-1], (1,len(b[0])-1))
    c_T = np.resize(c[0][0:len(c[0])-1], (1,len(c[0])-1))
    c_T[0][len(c_T[0])-1] = 0
    d_til = d.T[0:len(d[0])-1]
    v = np.zeros((1,len(a[0])-1))
    v[0][0] = a[0][0]
    v[0][len(v[0])-1] = c[0][len(c[0])-2]
    v = v.T
    
    if modo == 1:
        start_time = time.time()
        u,l = LU_decomp_tridiag(a_T,b_T,c_T)
        print("Vetor u =", u)
        print("vetor l =", l)
        print("Tempo decorrido para gerar solucao: ", time.time() - start_time)
        
        
    elif modo == 2:
        start_time = time.time()
        u,l = LU_decomp_tridiag(a_T,b_T,c_T)
        print("y_til =", resolve_sistema_nao_ciclico(u,l,c_T,d_til.T) )
        print("z_til= ", resolve_sistema_nao_ciclico(u,l,c_T,v.T))
        print("Tempo decorrido para gerar solucao: ", time.time() - start_time)
        
    elif modo == 3:
        start_time = time.time()
        print("x =", resolve_sistema_ciclico(a,b,c,d,n))
        print("Tempo decorrido para gerar solucao: ", time.time() - start_time)
   
            
        


  
    #r = np.dot(A,x) - d.T
    #print(r)
