# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 20:13:43 2018

@author: Michele97
"""


#test: A = [[1.0, -2.0, 3.0], [-1.0, 2.0, 1.0], [2.0, -1.0, -1.0]] 
#      b =[1.0, 0.0, -1.0]



from scipy import *
from numpy import *


def sisLin(A,b):
    [A, b] = gauss(A, b)
    n = len(b)
    
    #algoritmo di sostituzione all'indietro
    tol = 1e-15
    x = zeros(shape=(n,1))
    for i in range (n-1, -1, -1):
        ris = 0
        for j in range (i + 1, n):
            ris = ris + (A[i, j] * x[j])
        
        if abs(A[i,i]) < tol:
            print ("Errore!")
            return
        else:
            x[i] = (b[i] - ris)/ A[i, i]   
    
    return x
    



def gauss(A,b):
    [m, n] = shape(A)
    if(m != n):
        print ("errore! Matrice non quadrata!")
        return
    else:
        A = copy(A)
        b = copy(b)
        tol = 1e-15
        for i in range (0, n-1):
            if(A[i,i] == 0 or A[i,i] == tol):
                tecnicaPivot(A, b, i)
                            
            pivot = A[i,i]
                            
            #calcolo moltiplicatori
            for r in range (i+1, n):
                mol = (-1) * (A[r, i] / pivot)
                b[r] = (mol * b[i]) + b[r]

                for c in range (0, n):
                    A[r,c] = (mol * A[i,c]) + A[r,c]
                        
        return A, b
        
    
    
def tecnicaPivot(A, b, i):
    [m, n] = shape(A)
    tol = 1e-15
    trovato = False
    j = i + 1
    while (j in range (i+1, n) and not(trovato)):
        if(A[j,i] != tol and A[j,i] != 0):
            trovato = True
            for k in range (0, n):
                temp = A[i,k]
                A[i,k] = A[j,k]
                A[j,k] = temp
                            
            temp = b[i]
            b[i] = b[j]
            b[j] = temp
        j = j + 1
    
    