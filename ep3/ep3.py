import numpy as np
import math
import sys
sys.path.insert(0, '/home/fuoco/calculo_numerico/ep1')
from ep1 import resolve_sistema_nao_ciclico, LU_decomp_tridiag
sys.path.insert(0,'/home/fuoco/calculo_numerico/ep2')
from ep2 import integra, avalia_funcao, tabela_x, tabela_w

#Mateus Bonelli Salomão 11914789
#Pedro Pimentel Fuoco 11804313

def calcula_cs(n, h, k,x, w):                #funcao que calcula a diagonal superiopr
    c = []
    for i in range(1,n):
        c_i = (-1/(h**2))*integra(i*h, ((i+1)*h), k, x, w)
        c.append([c_i])
        #print('integral: ', integra(i*h, ((i+1)*h),  k, x, w))
        #print(c_i)
        #c.append(c_i)
    c.append([0])
    #return np.array(c)
    return np.array(c).T

def calcula_bs(n, h, k, x, w):                #calcula diagonal principal do sistema
    b = []
    for i in range(1,n+1):
        b_i = (1/(h**2))*((integra((i-1)*h, i*h, k, x, w))+integra(i*h, (i+1)*h, k, x, w))
        b.append([b_i])
        #b.append(b_i)
    return np.array(b).T
    #return np.array(b)

def calcula_as(n, h, k, x, w):                #funcao que calcula a diagonal inferior
    a = []
    a.append([0])
    for i in range(1,n):
        a_i = (-1/(h**2))*integra(i*h, ((i+1)*h), k, x, w) 
        a.append([a_i])
    return np.array(a).T


def calcula_ds(n, h, f, x, w):              #calcula lado direito do sistema
    d = []
    for i in range(1, n+1):
        d_i = (1/h)*(integra((i-1)*h, i*h, '(' +f + ')'+ '*(x-' + str((i-1)*h) +')', x, w) + integra(i*h, (i+1)*h, '('+ f +')'+'*(' +str((i+1)*h) +'-x)', x,w))      # construo a funcao de dentro concatenando os string e depois avalio com o avalia_vunao
        d.append([d_i])
        #d.append(d_i)
    return np.array(d).T
    #return np.array(d)

def calcula_alphas(n, h,f, k, x, w):
    a= calcula_as(n, h, k, x, w)
    b= calcula_bs(n, h, k, x, w)
    c= calcula_cs(n, h, k, x, w)
    d = calcula_ds(n, h, f, x, w)
    u, l = LU_decomp_tridiag(a,b,c)
    return resolve_sistema_nao_ciclico(u, l, c, d)
    
def calcula_erros(alphas, u, h, n):
    erros = []
    for i in range(0, n):
        #print(alphas[i])
        #print((i+1)*h)
        #print(avalia_funcao(u,(i+1)*h))
        erro_xi = abs(alphas[0][i]-avalia_funcao(u, (i+1)*h))
        erros.append(erro_xi)
    return np.array(erros)

def resolve_alphas(n, h,f, k, x, w):
    alpha = calcula_bs(n, h, k, x, w) #diag principal n elementos
    beta = calcula_cs(n, h, k, x, w)    #diags sup. e inferior - n-1 el,
    d = calcula_ds(n, h, f, x, w)       #lado direito n el.
    a = np.zeros(n)     #0 a n-1
    eta = np.zeros(n-1) #0 a n-2
    z = np.zeros(n) #0 a n-1
    a[0] = alpha[0]
    eta[0] = beta[0]/alpha[0]
    z[0] = d[0]/a[0]
    for i in range(1, n-2):
        a[i] = alpha[i]-beta[i-1]*eta[i-1]
        eta[i] = beta[i]/alpha[i]
        z[i] = (d[i]-beta[i-1]*z[i-1])/a[i]
    a[n-1] = alpha[n-1]-beta[n-2]*eta[n-2]
    z[n-1] = (d[n-1]-beta[n-2]*z[n-2])/a[n-1]
    alphas = np.zeros(n)
    alphas[n-1] = z[n-1]
    for i in range(n-2, -1, -1):
        alphas[i] = z[i] - eta[i]*alphas[i+1]
    return alphas

def calcula_erro_max(erros):
    return np.amax(erros)

def main():
    x1 = np.array([-math.sqrt(1/3), math.sqrt(1/3)])          #nos e pesos para a formula de Gauss com 2 pontos
    w1 = np.array([1, 1])
    n = int(input('Digite n: '))
    # intervalo [0,1]
    h = 1/(n+1)
    k = '1'
    f = '12*x*(1-x)-2'
    u = '(x**2)*((1-x)**2)'
    #k = '(math.e)**x'
    #f = '1 + (math.e)**(x)'
    #u = '(x-1)*(-1+(math.e)**(-x))'    

    alphas1 = calcula_alphas(n, h,f, k, x1, w1)
    erros1 = calcula_erros(alphas1, u, h, n)
    erro_max1 = calcula_erro_max(erros1)

    print("\nalpha_is1 encontrados: \n",alphas1)
    print("\nerros1 nos x_i: \n",erros1)
    print("\nerro1 maximo: \n",erro_max1)

if __name__=='__main__': 
    main()