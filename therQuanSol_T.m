function [therQuan] = therQuanSol_T(Z,T)
%{
    param Z: partition function array
    param T: temperature array
    return therQuan(Struct)
    {
        F: free energy          F = -kTlnZ
        U: internal energy      U = -{∂lnZ/∂β} = kT^2*{∂lnZ/∂T}
        Cv: specific heat       Cv = kβ^2{∂^2lnZ/∂β^2} = 2kT*{∂lnZ/∂T} + kT^2*{∂^2lnZ/∂T^2}
        S: entropy              S = k(lnZ - β*{∂lnZ/∂β}) = klnZ + kT*{∂lnZ/∂T}
    }

    k = 1 is boltzmann constant
%}
    %> Derivation
    lengthArr = length(T);
    logZ = log(Z);

    d_logZ_T = zeros(1,lengthArr);
    d_logZ_T(2:end-1) = (logZ(3:end)-logZ(1:end-2)) ./ (T(3:end)-T(1:end-2));

    d_logZ2_T2 = zeros(1,lengthArr);
    d_logZ2_T2(2:end-1) = (logZ(3:end)-2*logZ(2:end-1)+logZ(1:end-2)) ./ (T(3:end)-T(2:end-1)) ./ (T(2:end-1)-T(1:end-2));

    %> Free Energy: F
    F = - T .* log(Z);
    %> Internal Energy: U
    U = T.^2 .* d_logZ_T;
    U(1) = U(2);U(end) = U(end-1);
    %> Specific Heat: Cv
    Cv = 2 * T .* d_logZ_T + T.^2 .*d_logZ2_T2;
    Cv(1) = Cv(2);Cv(end) = Cv(end-1);
    %> Entropy: S
    S = logZ + T .* d_logZ_T;
    S(1) = S(2);S(end) = S(end-1);

    %> Save data
    therQuan = struct('F',F,'U',U,'Cv',Cv,'S',S);