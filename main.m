File_list = load('C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data/J2_0.5J1_h/File_list.mat').File_list;
File_len = length(File_list);

D = 20;
% normError = 1e-9;
normError = 5e-8;
tic;
for iT = 1:File_len

    loadpath = ['C:/Users/JAH/Desktop/组/张量网络/双层XY模型(暑期重做)/data/J2_0.5J1_h/',strtrim(File_list(iT,:)),'/'];
    Oloc = load([loadpath,'Oloc.mat']).Oloc;
    savepath = loadpath;

    umps = struct;
    umps.H = Oloc;
    VUMPS_single(umps,D,normError,savepath);

end
toc