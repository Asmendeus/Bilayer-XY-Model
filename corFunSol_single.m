function [x,y] = corFunSol_single(umps,M0,Mx,r)
%{
    param umps(struct)
    {
        name H: local tensor of 2D infinity TN
        name FL,FR: left-right environment tensors
        name AL,AR: left-right orthogonality tensors of the UMPS
        name AC,C: centeral tensors of the UMPS
        name Z: partition function of single site
    }
    param M0,Mx: corresponding impurity tensors at site = 0,x
    param r: maximum correlation length
    return x: correlation length array
    return y: correlation function value array
%}
%{
    With VUMPS algorithm, correlation function G(x) can be expressed as below:

             |‾‾|—— AC  —— AR  —— ··· —— AR  —— AR  ——|‾‾|      /          |‾‾|—— C  ——|‾‾|
             |  |   |       |             |     |     |  |     /           |  |        |  |
      G(x) = |FL|—— M0  ——  H  —— ··· ——  H  —— Mx  ——|FR|    /    Z^(x+1) |FL|————————|FR|
             |  |   |       |             |     |     |  |   /             |  |        |  |
             |__|—— AC* —— AR* —— ··· —— AR* —— AR* ——|__|  /              |__|—— C* ——|__|

             |‾‾|—— AR  —— ··· —— AR  ——|‾‾|      /          |‾‾|—— C  ——|‾‾|
             |  |    |             |    |  |     /           |  |        |  |
      G(x) = |F0|——  H  —— ··· ——  H  ——|Fx|    /    Z^(x+1) |FL|————————|FR|
             |  |    |             |    |  |   /             |  |        |  |
             |__|—— AR* —— ··· —— AR* ——|__|  /              |__|—— C* ——|__|

    where (x-1) local tensors H are between the impurity tensors M0 and Mx, and

      |‾‾|——   |‾‾|—— AC  ——     ——|‾‾|   —— AR  ——|‾‾|
      |  |     |  |   |            |  |      |     |  |
      |F0|—— = |FL|—— M0  ——  ,  ——|Fx| = —— Mx  ——|FR|
      |  |     |  |   |            |  |      |     |  |
      |__|——   |__|—— AC* ——     ——|__|   —— AR* ——|__|

%}
%> 1.Initialization
x = 1:1:r;
y = zeros(1,r);

%> 2.Solve correlation function y = G(x)
%· 1.generate F0,Fx
F0 = contractL_singleLayer(umps.FL,umps.AC,M0);
Fx = contractR_singleLayer(umps.FR,umps.AR,Mx);

%· 2.generate y = G(x)
y(1) = contract(F0,[1,2,3],Fx,[1,2,3]) / umps.Z^2;
(contract(contractL_singleLayer(F0,umps.AR,umps.H),[1,2,3],Fx,[1,2,3]) / contract(F0,[1,2,3],Fx,[1,2,3]))/umps.Z
for i = 2:1:r
    F0 = contractL_singleLayer(F0,umps.AR,umps.H) / umps.Z;
    y(i) = contract(F0,[1,2,3],Fx,[1,2,3]) / umps.Z^2;
end

tempDivisor = contract(umps.FL,[1],umps.C,[1]);
tempDivisor = contract(tempDivisor,[2],conj(umps.C),[1],[2,1,3]);
Divisor = contract(tempDivisor,[1,2,3],umps.FR,[1,2,3]);

y = y / Divisor;

    function Vout = contractL_singleLayer(Vin,AL,H)
        Vout = contract(Vin,[1],AL,[1]);
        Vout = contract(Vout,[1,3],H,[1,2]);
        Vout = contract(Vout,[1,4],conj(AL),[1,2]);
    end
    function Vout = contractR_singleLayer(Vin,AR,H)
        Vout = contract(AR,[3],Vin,[1]);
        Vout = contract(Vout,[2,3],H,[2,3]);
        Vout = contract(Vout,[4,2],conj(AR),[2,3]);
    end

end