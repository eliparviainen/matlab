
% use suitesparse
addpath /proj/matlab/local/toolbox/suitesparse
addpath /proj/matlab/local/toolbox/suitesparse/AMD/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/BTF/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/CAMD/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/CCOLAMD/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/CHOLMOD/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/COLAMD/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/CXSparse/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/KLU/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/LDL/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/LINFACTOR
addpath /proj/matlab/local/toolbox/suitesparse/MESHND
addpath /proj/matlab/local/toolbox/suitesparse/SPQR/MATLAB
addpath /proj/matlab/local/toolbox/suitesparse/SSMULT
addpath /proj/matlab/local/toolbox/suitesparse/UMFPACK/MATLAB


load examplegraph_nips
spv = ones(size(spr)); % unweighted

datestr(now, 'HH:MM:SS')

numG = 20;
walk_stop_thr = 1/20;
local_nei_thr = 0.01*walk_stop_thr;
numiter = 200;
gsneopt = defaultopts_gsne(numG, walk_stop_thr, local_nei_thr, numiter);
X = gsne(N, spr, spc, spv, gsneopt);
figure, blueyellowplot(X, spr, spc);
title(sprintf('GSNE example (V=%d, E=%d)',N,length(spr)),'fontsize',12);

datestr(now, 'HH:MM:SS')
