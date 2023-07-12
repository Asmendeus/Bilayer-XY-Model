function [C,Emax] = maxC_EffecHamilton(FL,FR)
%{
    param FL,FR: left-right environment tensor
    return C: central tensor of uMPS
%}
%{
    |‾‾|—— C ——|‾‾|
    |  |       |  |
    |FL|———————|FR| = Emax  —— C ——
    |  |       |  |
    |__|——   ——|__|

    here, Emax should be norm largest
%}

%> read shape data
FL_shape = size(FL);
D = FL_shape(1);

%> setting parameters of eigs
n = D*D;    % Heff will reshape into a n*n matrix
m = 1;      % number of maximum eigenvalues

opts.issym = false;
opts.v0 = rand(n,1);
opts.tol = 1e-14;

sigma = 'lm';    % largest magnitude

%> solve central tensor C
[maxVec, maxVal] = eigs(@CProduct,n,m,sigma,opts);
C = reshape(maxVec,[D,D]);
Emax = maxVal;

    %> subfunction: update C
    function Vout = CProduct(Vin)

        Vin = reshape(Vin,[D,D]);

        Vout = contract(FL,[1],Vin,[1]);
        Vout = contract(Vout,[3,1],FR,[1,2]);

        Vout = reshape(Vout,[n,1]);

    end

end