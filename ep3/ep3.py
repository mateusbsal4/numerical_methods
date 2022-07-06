import numpy as np
import math
import sys
sys.path.insert(0,'../ep1')
from ep1 import resolve_sistema_ciclico
sys.path.insert(0,'../ep2')
from ep2 import integra, avalia_funcao

#Mateus Bonelli Salomão 11914789
#Pedro Pimentel Fuoco 11804313


def calcula_as(n, h, k, x, w):                #funcao que calcula a diagonal inferior
    a = []
    a.append([0])
    for i in range(1,n):
        a_i = (-1/(h**2))*integra(i*h, ((i+1)*h), k, x, w) 
        a.append([a_i])
    return np.array(a).T


def calcula_cs(n, h, k,x, w):                #funcao que calcula a diagonal superiopr
    c = []
    for i in range(1,n):
        c_i = (-1/(h**2))*integra(i*h, ((i+1)*h), k, x, w) 
        c.append([c_i])
    c.append([0])
    return np.array(c).T


def calcula_bs(n, h, k, x, w):                #calcula diagonal principal do sistema
    b = []
    for i in range(1,n+1):
        b_i = (1/(h**2))*(integra((i-1)*h, i*h, k, x, w)+integra(i*h, (i+1)*h, k, x, w))
        b.append([b_i])
    return np.array(b).T


def calcula_ds(n, h, f, x, w):              #calcula lado direito do sistema
    d = []
    for i in range(1, n+1):
        d_i = (1/h)*(integra((i-1)*h, i*h, f + '*(x-' + str((i-1)*h) +')', x, w) + integra(i*h, (i+1)*h, str((i+1)*h) +'-x', x,w))      # construo a funcao de dentro concatenando os string e depois avalio com o avalia_vunao
        d.append([d_i])
    return np.array(d).T


def main():
    x = np.array([-math.sqrt(1/3), math.sqrt(1/3)])          #nos e pesos para a formula de Gauss com 2 pontos
    w = np.array([1, 1])
    n = int(input('Digite n: '))
    # intervalo [0,1]
    h = 1/(n+1)
    k = '1'
    f = '1'
    a= calcula_as(n, h, k, x, w)
    b= calcula_bs(n, h, k, x, w)
    c= calcula_cs(n, h, k, x, w)
    d  = calcula_ds(n, h, f, x, w)
    print("a: \n",a)
    print("b: \n",b)
    print("c: \n",c)
    print("d: \n",d)
    print("\nResolucao do sistema linear: \n", resolve_sistema_ciclico(a,b,c,d,n))


if __name__=='__main__': 
    main()
