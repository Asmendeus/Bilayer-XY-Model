import numpy as np
from scipy.special import iv
import os
import scipy.io as scio

#> 0.设置参量
N = 3       # Bessel函数取到的阶数
Nm = 2*N-1  # 8阶张量bond大小, M1 = M2 = M3 = M4 = N1 = N2 = N3 = N4 = 2*N-1
Ns = Nm**2  # reshape后4阶张量bond大小
J1 = 0.2         # 上层层内耦合强度
J2 = 0.5 * J1    # 下层层内耦合强度
K = 0.5 * J1     # 层间二阶Josephson耦合强度

T_list = np.linspace(0.5,1.3,81) * J1    # 温度
File_list = []

for T in T_list:

    beta = 1/T  # 倒温度, beta = 1/kT, 玻尔兹曼常数 k记为1

    #> 1.生成不等价张量Oloc
    #· 1.生成 m,n 的阶数向量
    order = np.zeros(Nm)
    for ii in range(N-1):
        order[ii*2+1] = ii+1
        order[ii*2+2] = -ii-1
    order = np.array(order,int)

    #· 2.生成根号下八项相乘对应生成的8阶张量
    IJ1 = iv(order,beta*J1)         # In(beta*J1)
    IJ2 = iv(order,beta*J2)         # Im(beta*J2)
    IJ1sq = np.sqrt(IJ1)
    IJ2sq = np.sqrt(IJ2)
    Isql = np.outer(IJ1sq,IJ2sq)
    Isq = np.outer(np.outer(np.outer(Isql,Isql),Isql),Isql).reshape(Nm,Nm,Nm,Nm,Nm,Nm,Nm,Nm) # 根号下的八项相乘对应生成的8阶张量

    #· 3.生成8阶张量的阶数索引张量
    indexm1 = np.array(np.outer(order,np.ones((Nm,Nm,Nm,Nm,Nm,Nm,Nm))),int).reshape(Nm,Nm,Nm,Nm,Nm,Nm,Nm,Nm)
    indexn1 = np.array(np.outer(np.ones((1,Nm)),np.outer(order,np.ones((Nm,Nm,Nm,Nm,Nm,Nm)))),int).reshape(Nm,Nm,Nm,Nm,Nm,Nm,Nm,Nm)
    indexm2 = np.array(np.outer(np.ones((Nm,Nm)),np.outer(order,np.ones((Nm,Nm,Nm,Nm,Nm)))),int).reshape(Nm,Nm,Nm,Nm,Nm,Nm,Nm,Nm)
    indexn2 = np.array(np.outer(np.ones((Nm,Nm,Nm)),np.outer(order,np.ones((Nm,Nm,Nm,Nm)))),int).reshape(Nm,Nm,Nm,Nm,Nm,Nm,Nm,Nm)
    indexm3 = np.array(np.outer(np.ones((Nm,Nm,Nm,Nm)),np.outer(order,np.ones((Nm,Nm,Nm)))),int).reshape(Nm,Nm,Nm,Nm,Nm,Nm,Nm,Nm)
    indexn3 = np.array(np.outer(np.ones((Nm,Nm,Nm,Nm,Nm)),np.outer(order,np.ones((Nm,Nm)))),int).reshape(Nm,Nm,Nm,Nm,Nm,Nm,Nm,Nm)
    indexm4 = np.array(np.outer(np.ones((Nm,Nm,Nm,Nm,Nm,Nm)),np.outer(order,np.ones((1,Nm)))),int).reshape(Nm,Nm,Nm,Nm,Nm,Nm,Nm,Nm)
    indexn4 = np.array(np.outer(np.ones((Nm,Nm,Nm,Nm,Nm,Nm,Nm)),order),int).reshape(Nm,Nm,Nm,Nm,Nm,Nm,Nm,Nm)

    #· 4_1.生成根号外部分对应的8阶张量
    delta1 = indexn3 + indexn4 - indexn1 - indexn2      # delta{n3+n4}{n1+n2+2k}对应的8阶索引张量,元素是2k在该位置的值
    delta2 = indexm1 + indexm2 - indexm3 - indexm4      # delta{m1+m2}{m3+m4+2k}对应的8阶索引张量,元素是2k在该位置的值

    deltak = delta1 - delta2
    deltak = np.array(deltak,dtype=bool)
    deltak = np.array(~deltak,dtype=int)              # delta1和delta2的2k值相等的位置置1, 不相等的位置置0
    deltakk = np.array(np.mod(delta1,2),dtype=bool)
    deltakk = np.array(~deltakk,dtype=int)            # k为整数的位置置1, 为半整数的位置置0

    Ik = iv(delta1/2,beta*K)        # Ik(beta*K)
    deltaIK = Ik*deltak*deltakk     #根号外的部分对应的8阶张量

    #· 5_1.生成不等价张量Oloc
    Oloc = Isq * deltaIK
    Oloc = Oloc.reshape(Ns,Ns,Ns,Ns)

    #· 4_1.生成M0,Mx根号外部分对应的8阶张量
    M0delta1 = indexn3 + indexn4 - indexn1 - indexn2 - 1    # delta{n3+n4}{n1+n2+2k+1}对应的8阶索引张量,元素是2k在该位置的值
    M0delta2 = indexm1 + indexm2 - indexm3 - indexm4        # delta{m1+m2}{m3+m4+2k}对应的8阶索引张量,元素是2k在该位置的值
    M0deltak = M0delta1 - M0delta2
    M0deltak = np.array(M0deltak,dtype=bool)
    M0deltak = np.array(~M0deltak,dtype=int)                # M0delta1和M0delta2的2k值相等的位置置1, 不相等的位置置0
    M0deltakk = np.array(np.mod(M0delta1,2),dtype=bool)
    M0deltakk = np.array(~M0deltakk,dtype=int)              # k为整数的位置置1, 为半整数的位置置0
    Ik = iv(M0delta1/2,beta*K)                              # Ik(beta*K)
    M0deltaIK = Ik*M0deltak*M0deltakk                       #根号外的部分对应的8阶张量

    Mxdelta1 = indexn3 + indexn4 - indexn1 - indexn2 - 1    # delta{n3+n4}{n1+n2+2k+1}对应的8阶索引张量,元素是2k在该位置的值
    Mxdelta2 = indexm1 + indexm2 - indexm3 - indexm4        # delta{m1+m2}{m3+m4+2k}对应的8阶索引张量,元素是2k在该位置的值
    Mxdeltak = Mxdelta1 - Mxdelta2
    Mxdeltak = np.array(Mxdeltak,dtype=bool)
    Mxdeltak = np.array(~Mxdeltak,dtype=int)                # Mxdelta1和Mxdelta2的2k值相等的位置置1, 不相等的位置置0
    Mxdeltakk = np.array(np.mod(Mxdelta1,2),dtype=bool)
    Mxdeltakk = np.array(~Mxdeltakk,dtype=int)              # k为整数的位置置1, 为半整数的位置置0
    Ik = iv(Mxdelta1/2,beta*K)                              # Ik(beta*K)
    MxdeltaIK = Ik*Mxdeltak*Mxdeltakk                       #根号外的部分对应的8阶张量

    #· 5_1.生成不等价张量M0,Mx
    M0 = Isq * M0deltaIK
    M0 = M0.reshape(Ns,Ns,Ns,Ns)

    Mx = Isq * MxdeltaIK
    Mx = Mx.reshape(Ns,Ns,Ns,Ns)

    #· 6.保存张量
    savepath = 'C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data/J2_0.5J1_h/T_' + str(round(T/J1,2)) + 'J1'
    os.mkdir(savepath)

    savepathOloc = savepath + '/Oloc.mat'
    savepathM0 = savepath + '/M0.mat'
    savepathMx = savepath + '/Mx.mat'

    savedictOloc = {'Oloc':Oloc}
    savedictM0 = {'M0':M0}
    savedictMx = {'Mx':Mx}

    scio.savemat(savepathOloc,savedictOloc)
    scio.savemat(savepathM0,savedictM0)
    scio.savemat(savepathMx,savedictMx)

    File_list.append('T_' + str(round(T/J1,2)) + 'J1')

savepathFile_list = 'C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data/J2_0.5J1_h/File_list.mat'
savepathT_list = 'C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data/J2_0.5J1_h/T_list.mat'
savedictFile_list = {'File_list':File_list}
savedictT_list = {'T_list':T_list}
scio.savemat(savepathFile_list,savedictFile_list)
scio.savemat(savepathT_list,savedictT_list)