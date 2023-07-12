% loadpath = 'C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data/J2_0.5J1_h/';
% File_list = load([loadpath,'File_list.mat']).File_list;
% File_len = length(File_list);
SE = zeros(1,81);
Z = zeros(1,81);

% T = load([loadpath,'T_list.mat']).T_list;

% for iT = 1:File_len
%     umps = load([loadpath,strtrim(File_list(iT,:)),'/umps.mat']).umps;
%     SE(iT) = entaEntrSol(umps.C);
%     Z(iT) = umps.Z;
% end
% [therQuan] = therQuanSol_T(Z,T);
% T = T / 0.2;

T = linspace(0.5,1.30,81);
for iT = 1:81
    umps = load(['C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data_SC/data/',num2str(iT+49),'/umps.mat']).umps;
    SE(iT) = entaEntrSol(umps.C);
    Z(iT) = umps.Z;
end
[therQuan] = therQuanSol_T(Z,T);

% plot(T,SE,'r',T,therQuan.F,'b',T,therQuan.U,'c',T,therQuan.Cv,'g',T,therQuan.S,'k');
% legend('SE','Free Energy','Internal Energy','Specific Heat','Entropy')
% plot(T,SE);
% title('SE~T');
% plot(T,therQuan.F);
% title('FreeEnergy~T');
% plot(T,therQuan.U);
% title('InternalEnergy~T');
% plot(T,therQuan.Cv);
% title('Cv~T');
% plot(T,therQuan.S);
% title('Entropy~T');

yyaxis left
plot(T,SE,'r');
ylabel('SE')
yyaxis right
plot(T,therQuan.Cv,'b');
ylabel('Cv')
