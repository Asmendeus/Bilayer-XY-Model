function AL = updateAL_LeftPolar(AC,C)
%{
    —— AC ——  left polar decompition  —— UL_AC —— PL_AC ——
       |      ======================       |        |

              left polar decompition
    —— C ——   ======================  —— UL_C —— PL_C ——

    —— AL ——  ===  —— UL_AC —— UL_C^{dag} ——
       |                |
%}

AC_shape = size(AC);
D = AC_shape(1);
d2 = AC_shape(2);

tempACL = reshape(AC,[D*d2,D]);
[U,~,V] = svd(tempACL,'econ');
UL_AC = reshape(U * V',[D,d2,D]);

[U,~,V] = svd(C,'econ');
UL_C = U*V';

AL = contract(UL_AC,[3],UL_C',[1]);