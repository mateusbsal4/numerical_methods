import numpy as np
import sympy as sym
import math
import sys
import time





#def acha_Legendren(n): #funcao que acha as raizes do polinomio de Legendre de ordem n+1, para n >= 2
#    x = sym.Symbol('x')
#    P_0 = 1
#    P_1 = x
#    Legendre = [P_0,P_1]
#    for i in range(2,n+1):
#        P_n = ((2*(i-1)+1)*x*Legendre[i-1] - (i-1)*Legendre[i-2])/(i)      
#        P_n = sym.simplify(P_n)
#        #print(P_n)
#        Legendre.append(P_n)
#        #print(Legendre)
#    return P_n
#
#def raizes_Legendren(Pn): #retorna uma tupla com os valores numéricos das raízes 
#    x = sym.Symbol('x')
#    return sym.solve(Pn,x)

#def acha_pesos(n,f):
#   W  = []
#   for i in range(n):
#       new_s = sym.Symbol("w"+str(i))
#       Symbols.append(new_s)
#   equations = []
#   for i in range(n+1):
#        equation = 0
#       for i in range(n):
#           equation += W[i]*f(raizes_Legendre(acha_Legendren(n))[i])
#        equation -= inacha_Legendren(i)
#    equations.append(equation) 

def tabela_x(n, L):                                         #devolve o array completo das raizes xi, dado n e as raizes positivas
    if n % 2 == 0:
        a = - L
        return np.concatenate((a, [0], L))                  #devolve o array completo das raizes xi, dado n e as raizes positivas
    return np.concatenate((a, L))

def tabela_w(n, L):                                         #cria a lista completa dos pesos wi dado n e os pesos positivos, 
    if n % 2 ==0:
        return np.concatenate((L, [2-2*(L.sum())], L))      #usa o fato de que a soma de todos os pesos é o tamanho do intervalo
    return np.concatenate((L, L))


def integra(a,b,f,n,x,w):
    assert len(x) == len(w)
    assert type(a) == float and type(b) == float
    I = 0
    if a == -1 and b == 1:                             #acha integral para o caso [a,b]= [1,1]
        for i in range(len(x)):
            I += w[i]*f(x[i])
        return I
    for i in range(len(x)):
        I += w[i]*f(((b-a)*x[i]+(b+a))/2)            #para um intervalo generico faz a mudança y=(2x-a-b)/(b-a) calcula a integral com y em [1,1]
    I *= (b-a)/2
    return I
    
    
if __name__ == "__main__":
    #dados do enunciado
    x6 = np.array([0.2386191860831969086305017, 0.6612093864662645136613996, 0.9324695142031520278123016])
    x8 = np.array([0.1834346424956498049394761, 0.5255324099163289858177390, 0.7966664774136267395915539, 0.9602898564975362316835609])	
    x10 = np.array([0.1488743389816312108848260, 0.4333953941292471907992659, 0.6794095682990244062343274, 0.8650633666889845107320967, 0.9739065285171717200779640])	
    w6 = np.array([0.4679139345726910473898703, 0.3607615730481386075698335, 0.1713244923791703450402961])
    w8 = np.array([0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314])
    w10 = np.array([0.2955242247147528701738930, 0.2692667193099963550912269, 0.2190863625159820439955349, 0.1494513491505805931457763, 0.0666713443086881375935688])

    start_time = time.time()
    print(integra(-1,1, lambda x: x**2, 6, tabela_x(6,x6), tabela_w(6,w6)))           #funcoes f de entrada devem ser do tipo lambda
    print(integra(-1,1, lambda x: (math.e)**x**2, 6, tabela_x(6,x6), tabela_w(6,w6)))
    print(integra(-1,1, lambda x: math.sin(x), 6, tabela_x(6,x6), tabela_w(6,w6)))
    print(integra(1,2, lambda x: 1/x, 6, tabela_x(6,x6), tabela_w(6,w6)))
    print("Tempo decorrido para gerar solucao: ", time.time() - start_time)


