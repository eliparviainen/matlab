function blueyellowplot(X, spr, spc, bw, nodots, msize)
%  
%  blueyellowplot(X, spr, spc, bw, nodots)
% 
%  Simple graph plotter.
%
%     X         Nx2 matrix, 2D-coordinates of nodes
%     spr, spc  graph links in sparse matrix format 
%     bw        set to 1 for grayscale (default 0, jet colormap)
%     nodots    set to 1 so node markers are not plotted (default 0)
%     msize     marker size (default 4)
%
%
% (c) Eli Parviainen, 2014
% Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
% 
  
  
  % default is yellowish to blue plus red markers, bw gives gray to black
  if ~exist('bw','var') | isempty(bw), bw=0; end;  
  
  % default is to plot markers at nodes
  if ~exist('nodots','var') | isempty(nodots), nodots=0; end;
  
  if ~exist('msize','var') | isempty(msize), msize=4; end;  
  
  clf
  hold on

  linkpairs = [spr(:) spc(:)];
  numcolors = 8;
  
  if bw>1, numcolors=1, end;
  
  N=size(X,1);
  
  hold on

  % scale to 0,100 so color settings work at least somehow
  scX = bsxfun(@minus, X, min(X));
  scX = bsxfun(@rdivide, scX, max(scX)); 
  scX = 100*scX;

  % this matters when printing on paper
  linewid=0.1;
  
  
  
  %linkpairs = linkpairs(find(linkpairs(:,1)>0&linkpairs(:,2)>0),:);

  % order from longest to shortest
  linkpairlen2 = sum((scX(linkpairs(:,1),:)-scX(linkpairs(:,2),:)).^2,2);
  [linkpairlen2 plotjarj] = sort(linkpairlen2,'descend');
  linkpairs = linkpairs(plotjarj,:);


  % logarith makes colors more readable
  cval=double(log(linkpairlen2+1));
  
  % a color index for each link pair
  cval=cval-min(cval);
  cval=cval/max(cval);
  cval=floor(cval*(numcolors-1))+1;
  
  if bw
    colors = gray(numcolors+2);
  else
    colors = jet(numcolors+3);
  end;

  % plot from longest to shortest
  for ind = 1:length(linkpairs)    
    line([X(linkpairs(ind,1),1), X(linkpairs(ind,2), 1)], ...
         [X(linkpairs(ind,1),2), X(linkpairs(ind,2), 2)], ...
         'color',colors(cval(ind),:), ...
         'linewidth',linewid ...
         );
  end;
  
  % dots at nodes
  if ~nodots
    if bw
      plot(X(:,1),X(:,2),'k.','markersize',msize);
    else
      plot(X(:,1),X(:,2),'r.','markersize',msize);
    end;
  end;
  
  axis tight
  axis off
  box off
  