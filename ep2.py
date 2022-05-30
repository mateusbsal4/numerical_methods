import numpy as np
import sympy as sym
import math
import sys
import time





def acha_Legendren(n): #funcao que acha as raizes do polinomio de Legendre de ordem n+1, para n >= 2
    x = sym.Symbol('x')
    P_0 = 1
    P_1 = x
    Legendre = [P_0,P_1]
    for i in range(2,n+1):
        P_n = ((2*(i-1)+1)*x*Legendre[i-1] - (i-1)*Legendre[i-2])/(i)      
        P_n = sym.simplify(P_n)
        #print(P_n)
        Legendre.append(P_n)
        #print(Legendre)
    return P_n

def raizes_Legendren(Pn): #retorna uma tupla com os valores numéricos das raízes 
    x = sym.Symbol('x')
    return sym.solve(Pn,x)


if __name__ == "__main__":
    print(acha_Legendren(2))
    print(acha_Legendren(5))
    print(raizes_Legendren(acha_Legendren(2)))
    


