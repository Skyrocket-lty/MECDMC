# -*- codeing = utf-8 -*-
# @Time : 2023/3/24 21:17
# @Author : 刘体耀
# @File : CMCFGRBF.py
# @Software: PyCharm


import numpy as np
#import pandas as pd
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import copy

# SM = np.loadtxt(r"SM_GRB.txt", dtype=float)
# MM = np.loadtxt(r"miRNA_GRB.txt", dtype=float)

SM = np.loadtxt(r"SM相似性矩阵.txt", dtype=float)
MM = np.loadtxt(r"miRNA相似性矩阵.txt", dtype=float)
Y1= np.loadtxt(r"SM-miRNA关联矩阵.txt",dtype=float)
SM_miRNA_k = np.loadtxt(r"SM-miRNA-已知关联.txt",dtype=int)
SM_miRNA_uk = np.loadtxt(r"SM-miRNA-未知关联.txt",dtype=int)




def SVT(Y, b, lty):
    #奇异值分解
    U,S,V= np.linalg.svd(Y)
    #print(S)
    #奇异值收缩
    for index in range(0,S.size):#s.size函数是矩阵元素的个数
        if S[index] >= b*lty[index]:
            S[index] = S[index] - b*lty[index]
        else:
            S[index] = 0


    #奇异值矩阵标准化
    s = np.diag(S)
    row , col = Y.shape[0] , Y.shape[1]
    if row < col:
        s_n = np.column_stack((s, np.zeros((row, col - row))))
    else:
        s_n = np.row_stack((s, np.zeros((row-col, col))))
    #U*奇异值收缩矩阵*V
    Y_n = np.dot(U, np.dot(s_n, V))
    return Y_n


def TNNR(alpha, beta,gamma, t, omega, tol1, tol2, maxiter,A,B,p,r):
    X = t
    W = X
    Y = X
    iter0 = 1
    stop1 = 1
    stop2 = 1

    # the processing of computing Wt
    M=np.dot(B.transpose(),A)
    U,S,V= np.linalg.svd(M)#这里面奇异值的个数只有一个
    u,s,v= np.linalg.svd(X)
    Wt=np.zeros_like(s)
    for i in range(0,S.size):
        Wt[i]=p*(1-S[i])*pow(s[i],p-1)
    if(Wt[i]<0):
        Wt[i]=0

    while stop1 > tol1 or stop2 > tol2:



        # the processing of computing W
        tran = (1/beta) * (Y+alpha*(t*omega))+X
        W = tran - (alpha/(alpha+beta))*omega*tran
        W[W < 0] = 0
        W[W > 1] = 1

        # the processing of computing X
        X_1 = SVT(W-(1/beta)*Y, 1/beta,Wt)

        # the processing of computing Y
        Y = Y + gamma*(X_1-W)

        stop1_0 = stop1
        if np.linalg.norm(X) != 0:
            stop1 = np.linalg.norm(X_1-X) / np.linalg.norm(X)
        else:
            stop1 = np.linalg.norm(X_1-X)
        stop2 = np.abs(stop1-stop1_0)/(max(1, np.abs(stop1_0)))
        X = X_1

        if iter0 >= maxiter:
            iter0 = maxiter
            print('reach maximum iteration,did not converge!')
            break
        iter0 = iter0 + 1
    T_recover = W
    return T_recover, iter0


def run_MC(t):
    # BNNR parameter
    maxiter = 500
    alpha = 2
    beta = 10
    gamma = 1
    p=0.8
    tol1 = 2 * 1e-3
    tol2 = 1 * 1e-5
    omega = np.zeros(t.shape)
    omega[t.nonzero()] = 1
    #插入第一层循环，或者叫第一步
    for i in range(0,2):
        U, S, V = np.linalg.svd(t)
        r =5
        A = U[:r, :]
        #print(np.dot(A,A.transpose()))

        B = V[:r, :]
        #print(B.shape)
        #print(np.dot(B,B.transpose()))

        t, k = TNNR(alpha, beta,gamma, t, omega, tol1, tol2, maxiter,A,B,p,r)
    Smmi = t
    return Smmi


#python


def TCMF(alpha, beta,gamma, Y, maxiter,A,B,C,SM,MM):

    iter0=1
    while True:

        a = np.dot(Y,B)+beta*np.dot(SM,A)
        b = np.dot(np.transpose(B),B)+alpha*C+beta*np.dot(np.transpose(A),A)
        c = np.dot(np.transpose(Y),A)+gamma*np.dot(MM,B)
        d = np.dot(np.transpose(A), A) + alpha * C + gamma * np.dot(np.transpose(B), B)

        A = np.dot(a,np.linalg.inv(b))
        B = np.dot(c, np.linalg.inv(d))

        if iter0 >= maxiter:

            #print('reach maximum iteration!')
            break
        iter0 = iter0 + 1
    Y= np.dot(A,np.transpose(B))
    Y_recover = Y
    return Y_recover


def run_MC_2(Y):
    maxiter = 1000
    alpha = 0.25
    beta = 0.01
    gamma = 0.01
    #SVD

    U, S, V = np.linalg.svd(Y)
    S=np.sqrt(S)
    r = 6
    Wt = np.zeros([r,r])
    for i in range(0,r):
        Wt[i][i]=S[i]
    U= U[:, 0:r]
    V= V[0:r,:]
    A = np.dot(U,Wt)
    B1 = np.dot(Wt,V)
    B=np.transpose(B1)
    C = np.zeros([r,r])
    for i in range(0,r):
        C[i][i] = 1
    #print(C)

    Y = TCMF(alpha, beta,gamma,Y, maxiter,A,B,C,SM,MM)
    Smmi = Y
    return Smmi



if __name__ == "__main__":

    X_1 = run_MC(Y1)  #Truncated schattenp norm minimization
    M_1 = run_MC_2(X_1)  #Low-rank matrix factorization


