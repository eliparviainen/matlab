
load examplegraph_nips
spv = ones(size(spr)); % unweighted

datestr(now, 'HH:MM:SS')

numG = 20;
walk_stop_thr = 1/20;
local_nei_thr = 0.01*walk_stop_thr;
numiter = 200;
gsneopt = defaultopts_gsne(numG, walk_stop_thr, local_nei_thr, numiter);
X = gsne_oldsparse(N, spr, spc, spv, gsneopt);
figure, blueyellowplot(X, spr, spc); 
title(sprintf('GSNE example (V=%d, E=%d)',N,length(spr)),'fontsize',12);

datestr(now, 'HH:MM:SS')
