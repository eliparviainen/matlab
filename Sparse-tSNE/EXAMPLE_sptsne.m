
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


load exampledata_mnist2000
x = double(x)/255;
y = double(y);

K = 10;
L = 4*K;
G = 40;
numiter = 1000;
starttime = datestr(now, 'HH:MM:SS')
X = sptsne(x, K, L, G, numiter);
endtime = datestr(now, 'HH:MM:SS')

figure, hold on
colors = hsv(max(y));
for ind = 1:max(y)
  cpts = find(y==ind);      
  h = plot(X(cpts,1), X(cpts,2), ...
           'marker', 'o', ...
           'linestyle','none', ...
           'markerfacecolor', colors(ind, :), ...
           'markeredgecolor', colors(ind,:), ...
           'markersize', 5 ...
           );
end;
axis tight
set(gca,'xtick',[],'ytick',[]);

title(sprintf('sparse t-SNE example (N=%d, K=%d, L=%d, G=%d)',size(x,1),K,L,G),'fontsize',12);


