function [AC,Emax] = maxAC_EffecHamilton(FL,FR,H)
%{
    param FL,FR: left-right environment tensor
    param H: local tensor of infinite 2D TN
    return AC: central tensor of uMPS
%}
%{
    |‾‾|—— AC ——|‾‾|
    |  |   |    |  |
    |FL|—— H  ——|FR| = Emax  —— AC ——
    |  |   |    |  |            |
    |__|——    ——|__|

    here, Emax should be norm largest
%}

%> read shape data
FL_shape = size(FL);
H_shape = size(H);
D = FL_shape(1);
d2 = H_shape(2);

%> setting parameters of eigs
n = D*d2*D;     % Heff will reshape into a n*n matrix
m = 1;          % number of maximum eigenvalues

opts.issym = false;
opts.v0 = rand(n,1);
opts.tol = 1e-14;

sigma = 'lm';    % largest magnitude

%> solve central tensor AC
[maxVec, maxVal] = eigs(@ACProduct,n,m,sigma,opts);
AC = reshape(maxVec,[D,d2,D]);
Emax = maxVal;

    %> subfunction: update AC
    function Vout = ACProduct(Vin)

        Vin = reshape(Vin,[D,d2,D]);

        Vout = contract(FL,[1],Vin,[1]);
        Vout = contract(Vout,[1,3],H,[1,2]);
        Vout = contract(Vout,[2,3],FR,[1,2]);

        Vout = reshape(Vout,[n,1]);

    end

end