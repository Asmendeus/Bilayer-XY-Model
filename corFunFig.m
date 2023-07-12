umps = load('C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data/J2_0.5J1_h/T_0.6J1/umps.mat').umps;
M0 = load('C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data/J2_0.5J1_h/T_0.6J1/M0.mat').M0;
Mx = load('C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data/J2_0.5J1_h/T_0.6J1/Mx.mat').Mx;
[x,y] = corFunSol_single(umps,M0,Mx,50);
plot(x,y);