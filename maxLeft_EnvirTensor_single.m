function [FL,ZL] = maxLeft_EnvirTensor_single(AL,H)
%{
    param AL: left orthogonality tensor of uMPS
    param H: local tensor of infinite 2D TN
    return FL: left environment tensor
    return ZL: max eigenvalue
%}
%{
    |‾‾|—— AL  ——      |‾‾|——
    |  |    |          |  |
    |FL|——  H  —— = ZL |FL|——     update: AL -> H ->AL*
    |  |    |          |  |
    |__|—— AL* ——      |__|——
%}

%> read shape data
AL_shape = size(AL);
H_shape = size(H);
D = AL_shape(1);
d1 = H_shape(1);

%> setting parameters of eigs
n = D*d1*D;     % TL will reshape into a n*n matrix
m = 1;          % number of maximum eigenvalues

opts.issym = false;
opts.v0 = rand(n,1);
opts.tol = 1e-14;

sigma = 'lm';    % largest magnitude

%> solve maximum eigenvalue and left environment tensor FL
[maxVec, maxVal] = eigs(@UpdateFL_single,n,m,sigma,opts);
FL = reshape(maxVec,[D,d1,D]);
ZL = maxVal;

    %> subfunction: update FL
    function Vout = UpdateFL_single(Vin)

        Vout = reshape(Vin,[D,d1,D]);

        Vout = contract(Vout,[1],AL,[1]);
        Vout = contract(Vout,[1,3],H,[1,2]);
        Vout = contract(Vout,[1,4],conj(AL),[1,2]);

        Vout = reshape(Vout,[n,1]);

    end

end