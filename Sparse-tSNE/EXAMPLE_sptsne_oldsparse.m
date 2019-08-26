load exampledata_mnist2000
x = double(x)/255;
y = double(y);

K = 10;
L = 4*K;
G = 40;
numiter = 1000;
starttime = datestr(now, 'HH:MM:SS')
X = sptsne_oldsparse(x, K, L, G, numiter);
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


