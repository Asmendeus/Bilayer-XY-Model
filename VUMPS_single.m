function [umps] = VUMPS_single(umps,D,normError,savepath)
%{
    param umps(struct):
    {
        name H: local tensor of 2D infinity TN
        name FL,FR: left-right environment tensors
        name AL,AR: left-right orthogonality tensors of the UMPS
        name AC,C: centeral tensors of the UMPS
        name Z: partition function of single site
    }
    param D: dimension of local tensor's virtual index
    param normError: the maximum error in judging convergence
    param savepath: savepath
    return umps(struct): as above
%}
%{
    Introduction of VUMPS algorithm
    =======================================
    The 2D infinity tensor network as below:
            ···  ···  ···  ···
             |    |    |    |
      ··· —— H —— H —— H —— H —— ···
             |    |    |    |
      ··· —— H —— H —— H —— H —— ···
             |    |    |    |
      ··· —— H —— H —— H —— H —— ···
             |    |    |    |
      ··· —— H —— H —— H —— H —— ···
             |    |    |    |
            ···  ···  ···  ···
    we marked as "H" for the local tensor with four legs:
             (d2)                                    2
              |                                      |
      (d1) —— H —— (d1)     order of 4 indexes: 1 —— H —— 3
              |                                      |
             (d2)                                    4
    =======================================
    State tensor for a single site: 3-order tensor (2 virtual + 1 physical).
    Marked as "A" for "State":
            (d2)                                  2
             |                                    |
      (D) —— A —— (D)    order of 3 indexes: 1 —— A —— 3
    The same for its Hermitian conjugate:
      (D) —— A* —— (D)                       1 —— A* —— 3
             |                                    |
            (d2)          order of 3 indexes:     2
    =======================================
    The target is to solve the variational boundary UMPS, meeting

             |    |    |    |
      ··· —— H —— H —— H —— H —— ···
             |    |    |    |                    |    |    |    |
      ··· —— A —— A —— A —— A —— ···  ∝  ··· —— A —— A —— A —— A —— ···
    =======================================
    In VUMPS, we solve the equivalent regular representation of the UMPS:

              |    |    |    |
       ··· —— A —— A —— A —— A —— ···

              |     |     |     |
    =  ··· —— AL —— AL —— AC —— AR —— ···

              |     |          |
    =  ··· —— AL —— AL —— C —— AR —— ···

    where AL/AR satisfy the left/right orthogonal condition.
    =======================================
%}

    %> 0.Check the parameters
    H_ndims = ndims(umps.H);
    if (H_ndims ~= 4)
        error('局域张量不是4-legs型, 请检查后重试!');
    end
    H_shape = size(umps.H);      % matlab的size = numpy的shape; matlab的numel = numpy的size
    if (H_shape(1) ~= H_shape(3))
        error('局域张量左右指标维数不一致, 请检查后重试!')
    end
    if (H_shape(2) ~= H_shape(4))
        error('局域张量上下指标维数不一致, 请检查后重试!')
    end
    if (D <= 0)
        error('虚拟指标维数必须为正整数, 请检查后重试!')
    end
    if (normError <= 0)
        error('收敛误差必须为正实数, 请检查后重试!')
    end

    %> 1.Get physical index dimention d for local tensor T and set other parameters
    d2 = H_shape(2);    % uMPS的物理指标维数, 局域张量上下指标维数
    num_iter = 0;               % 迭代次数
    epsilon = 1 + normError;    % 迭代误差

    %> 2.Initialize left-right orthogonal tensors AL,AR
    umps.AL = rand(D,d2,D);
    umps.AL = umps.AL / norm(umps.AL(:));
    umps.AR = rand(D,d2,D);
    umps.AR = umps.AR / norm(umps.AR(:));

    while (epsilon > normError)

        %> 3.Generate environment tensors FL,FR
        [umps.FL,ZL] = maxLeft_EnvirTensor_single(umps.AL,umps.H);
        [umps.FR,ZR] = maxRight_EnvirTensor_single(umps.AR,umps.H);

        umps.Z = (ZL+ZR) / 2;    % 收敛时, 应有 ZL = ZR = Z

        %> 4.Generate central tensors AC,C based on effective Hamiltonian
        [umps.AC,~] = maxAC_EffecHamilton(umps.FL,umps.FR,umps.H);
        [umps.C,~] = maxC_EffecHamilton(umps.FL,umps.FR);

        %> 5.Update the AL,AR based on AC,C
        umps.AL = updateAL_LeftPolar(umps.AC,umps.C);
        umps.AR = updateAR_RightPolar(umps.AC,umps.C);

        %> 6.Calculate the error and decide whether to continue the iteration
        errorL = umps.AC - contract(umps.AL,[3],umps.C,[1]);
        errorR = umps.AC - contract(umps.C,[2],umps.AR,[1]);
        normL = norm(errorL(:));
        normR = norm(errorR(:));
        epsilon = max(normL,normR);

        num_iter = num_iter + 1;
        fprintf('当前已迭代 %d',num_iter);fprintf(' 次\n');
        fprintf('误差为 %.4e',epsilon);fprintf('\n\n')

        SC = svd(umps.C)';
        SE = sum(-SC.^2 .* log(SC.^2));
        fprintf('纠缠熵为 %.4e',SE);fprintf('\n');
        fprintf('奇异谱为\n')
        disp(SC);

        %> Save umps (and ZL,ZR)
        save([savepath,'umps'],'umps')
        % save([savepath,'ZL'],'ZL')
        % save([savepath,'ZR'],'ZR')

    end