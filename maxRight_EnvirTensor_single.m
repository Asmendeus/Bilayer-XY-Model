function [FR,ZR] = maxRight_EnvirTensor_single(AR,H)
%{
    param AR: right orthogonality tensor of uMPS
    param H: local tensor of infinite 2D TN
    return FR: right environment tensor
    return ZR: max eigenvalue
%}
%{
    —— AR  ——|‾‾|      ——|‾‾|
        |    |  |        |  |
    ——  H  ——|FR| = ZR ——|FR|
        |    |  |        |  |
    —— AR* ——|__|      ——|__|
%}

%> read shape data
AR_shape = size(AR);
H_shape = size(H);
D = AR_shape(1);
d1 = H_shape(1);

%> setting parameters of eigs
n = D*d1*D;     % TR will reshape into a n*n matrix
m = 1;          % number of maximum eigenvalues

opts.issym = false;
opts.v0 = rand(n,1);
opts.tol = 1e-14;

sigma = 'lm';    % largest magnitude

%> solve maximum eigenvalue and right environment tensor FR
[maxVec, maxVal] = eigs(@UpdateFR_single,n,m,sigma,opts);
FR = reshape(maxVec,[D,d1,D]);
ZR = maxVal;

    %> subfunction: update FR
    function Vout = UpdateFR_single(Vin)

        Vout = reshape(Vin,[D,d1,D]);

        Vout = contract(AR,[3],Vout,[1]);
        Vout = contract(Vout,[2,3],H,[2,3]);
        Vout = contract(Vout,[4,2],conj(AR),[2,3]);

        Vout = reshape(Vout,[n,1]);

    end

end