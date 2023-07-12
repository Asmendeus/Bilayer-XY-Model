function AR = updateAR_RightPolar(AC,C)
%{
    —— AC ——  right polar decompition  —— PR_AC —— UR_AC ——
       |      =======================                |

              right polar decompition
    —— C ——   =======================  —— PR_C —— UR_C ——

    —— AR ——  ===  —— UR_C^{dag} —— UR_AC ——
       |                              |
%}

AC_shape = size(AC);
D = AC_shape(1);
d2 = AC_shape(2);

tempACR = reshape(AC,[D,d2*D]);
[U,~,V] = svd(tempACR,'econ');
UR_AC = reshape(U * V',[D,d2,D]);

[U,~,V] = svd(C,'econ');
UR_C = U*V';

AR = contract(UR_C',[2],UR_AC,[1]);