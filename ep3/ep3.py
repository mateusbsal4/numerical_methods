import numpy as np
import math
#from ep2 import integra, avalia_funcao 

def calcula_as(n, h, k, x, w):                #funcao que calcula a diagonal inferior
    a = []
    for i in range(1,n):
        a_i = (-1/(h**2))*integra(i*h, ((i+1)*h), k, x, w) 
        a.append([a_i])
    return np.array(a)




def integra(a,b,f,x,w):
    assert len(x) == len(w) 
    #assert type(a) == float and type(b) == float
    I = 0
    if a == -1 and b == 1:                             #acha integral para o caso [a,b]= [1,1]
        for i in range(len(x)):
            I += w[i]*avalia_funcao(f,(x[i]))
        return I
    for i in range(len(x)):
        I += w[i]*avalia_funcao(f,((b-a)*x[i]+(b+a))/2)            #para um intervalo generico faz a mudança y=(2x-a-b)/(b-a) calcula a integral com y em [1,1]
    I *= (b-a)/2
    return I

def avalia_funcao(f, y):
    x = y
    return eval(f)


def calcula_cs(n, h, k,x, w):                #funcao que calcula a diagonal superiopr
    return calcula_as(n, h, k, x, w)


def calcula_bs(n, h, k, x, w):                #calcula diagonal principal do sistema
    b = []
    for i in range(1,n+1):
        b_i = (1/(h**2))*(integra((i-1)*h, i*h, k, x, w)+integra(i*h, (i+1)*h, k, x, w))
        b.append([b_i])
    return np.array(b)


def calcula_ds(n, h, f, x, w):              #calcula lado direito do sistema
    d = []
    for i in range(1, n+1):
        d_i = (1/h)*(integra((i-1)*h, i*h, f + '*(x-' + str((i-1)*h) +')', x, w) + integra(i*h, (i+1)*h, str((i+1)*h) +'-x', x,w))      # construo a funcao de dentro concatenando os string e depois avalio com o avalia_vunao
        d.append([d_i])
    return np.array(d)
























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
    print(a)
    print(b)
    print(c)
    print(d)
    

if __name__=='__main__': main()
