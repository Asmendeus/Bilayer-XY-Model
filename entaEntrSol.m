function SE = entaEntrSol(C)

    SC = svd(C);
    SE = sum(-SC.^2 .* log(SC.^2));

end