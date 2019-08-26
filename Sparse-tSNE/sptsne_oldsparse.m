function X = sptsne_oldsparse(x, K, L, G, maxiter, verbose, tsneopt, sigma2)

%  -----------------------------------------------------------
%  This is a version of sptsne.m that uses original 
%  Matlab sparse. The version with SuitSparse sparse2 function
%  is faster and needs less memory.
%  -----------------------------------------------------------
%
%  X = sptsne(x, K, L, G, maxiter, tsneopt, sigma2)
%
%  Runs sparse t-SNE on the given data set. 
% 
%  IN:
%  x                  data, NxP
%  K                  perplexity aka effective #neighbors
%                     (this has no effect if sigma2 is given)
%  L                  number of local links per point, e.g. L=4K
%  G                  number of global links per point
%  maxiter (optional) number of iterations, default 1000
%  verbose (optional) 0=silent, 1=report progress (default)
%  tsneopt (optional) optimizer internal parameters (see code)
%  sigma2 (optional)  Nx1, sigma squared, pre-determined kernel widths
%                     (default is to call apprsigmas.m to get these)
%
%  OUT:
%  X                  Nx2, 2D coordinates
%                     example usage: plot(X(:,1), X(:,2), '.');
%
%
% (c) Eli Parviainen, 2014
% Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
% 

  
% ------------------
% init
% ------------------
  
  if ~exist('tsneopt', 'var') | isempty(tsneopt) 
    verbose = 1;
  end;
  
  N=size(x,1);
  x = double(x);


  if ~exist('tsneopt', 'var') | isempty(tsneopt)        
    tsneopt=struct();

    tsneopt.verbose = verbose;
    tsneopt.verboseiter = 50;
    tsneopt.max_iter = 1000; % may change below

    tsneopt.useEarlyExaggeration = 1;
    tsneopt.stop_lying_iter = 200;
    tsneopt.exaggFactor = 4;  

    tsneopt.outputdims = 2;
    tsneopt.momentum = 0.5;                                     
    tsneopt.final_momentum = 0.8;                               
    tsneopt.mom_switch_iter = 250;                             
    tsneopt.epsilon = 500;                                     
    tsneopt.min_gain = .01;                                    

    tsneopt.plotprogress = 0; 

  end; % if ~exist tsneopt


  % separately since changed more often than other tsneopts
  if exist('maxiter', 'var') & ~isempty(maxiter)
    tsneopt.max_iter = maxiter;
  end;

  
  % ------------------
  % link structure
  % ------------------

  % link locations
  
  if verbose, fprintf('\nsp_tsne: L- and G-links'); end;
  [spr spc spw LneiD2] = lggraph(x, L, G);

  spr=uint32(spr);
  spc=uint32(spc);
  spw=uint8(spw); % 1=local link, 2=global link


  % (this is logically not part of "link structure" but
  % is done here so LneiD2 can be cleaned away)
  if ~exist('sigma2', 'var') | isempty(sigma2)
    if verbose, fprintf('\nsp_tsne: sigmas'); end;
    sigma2 = apprsigmas(LneiD2,K);  
  end; % sigmas
  clear LneiD2 

  
  % neighborhood probabilities (spvv) for the existing links
  
  if verbose, fprintf('\nsp_tsne: weights for links'); end;
  
  % for highdim x runs out of memory if spr/spc is long
  % and x gets indexed by them, that's why batches
  d2 = zeros(length(spr)+length(spc),1);
  batch=10000;
  offs = 0:batch:length(spr);
  for bind=1:length(offs)
    range = offs(bind)+1:min(length(spr), offs(bind)+batch);
    d2(range) = sum((double(x(spr(range),:)) - double(x(spc(range),: ))).^2,2);
  end;
  offs = 0:batch:length(spc);
  for bind=1:length(offs)
    range = offs(bind)+1:min(length(spc), offs(bind)+batch);
    d2(length(spr)+range) = sum((double(x(spc(range),:)) - double(x(spr(range),: ))).^2,2);
  end;

  spvv = exp(-d2./sigma2([spr; spc]));
  clear d2 batch offs range
  clear x sigma2

  % ------------------
  % well-form similarities
  % ------------------  
  
  % make matrix rows sum to one by dividing by
  % the normalization factor sparseZ
  sparseZ = zeros(N,1);
  for ind=1:length(spr)  
    sparseZ(spr(ind)) = sparseZ(spr(ind)) + spvv(ind);
  end;
  for ind=1:length(spc)  
    sparseZ(spc(ind)) = sparseZ(spc(ind)) + spvv(length(spr)+ind);
  end;

  spv_norm = spvv./sparseZ([spr; spc]);
  spv_norm(isnan(spv_norm))=eps;
  clear spvv sparseZ

  % make the similarities symmetric
  % (the mess below symmetrizes the matrix without ever
  % creating a full matrix, to save some memory as compared
  % to the straightforward symM = (M+M')/2 symmetrization).

  if verbose, fprintf('\nsp_tsne: make symmetric...'); end;

  % here: no assumptions about order of link lists spr,spc,spv,spw
  sprYA = [spr; spc];
  spcYA = [spc; spr];
  clear spr spc
  spwYA = [spw; spw];
  clear spw
  
  % divide indices and lists into upper (y) and lower (a) tri
  yind = find(spcYA>sprYA);
  aind = find(spcYA<sprYA);

  % upper tri is named spr since they are the final values
  spr = sprYA(yind); sprA = sprYA(aind); clear sprYA
  spc = spcYA(yind); spcA = spcYA(aind); clear spcYA
  spvY = spv_norm(yind); spvA = spv_norm(aind); clear spv_norm aind
  spw = spwYA(yind); clear spwYA clear yind

  spr=uint32(spr);
  spc=uint32(spc);
  
  % lower tri
  % MATLAB SPARSE
  A = sparse(double(sprA),double(spcA),spvA+1,N,N);
  clear sprA spcA spvA

  % read lower tri values in col1 col2 ... order but call a column "row"
  [spcA sprA spvAp1] = find(A); 
  clear A

  % put lower tri value to upper tri in order row1, row2, ...
  % MATLAB SPARSE
  A = sparse(double(sprA),double(spcA),spvAp1,N,N);
  % read former lower tri from the upper tri, now the values have upper tri order
  clear spcA sprA
  [~,~, spvAp1] = find(A); 
  clear A

  % upper tri and reordered lower tri can be summed
  spv = (spvY(:)+spvAp1(:)-1)/2; clear spvY spvAp1 
  
  % here: sp=sparse([spr; spc], [spc; spr], [spv; spv], N, N);
  % would create a symmetric matrix

  % ------------------
  % dimred
  % ------------------  

  % separate local and global links

  lokind_flat = find(spw==1);
  globind_flat = find(spw==2);
  clear spw

  sprL = spr(lokind_flat);
  spcL = spc(lokind_flat);
  spvL = spv(lokind_flat);
  clear lokind_flat

  sprG = spr(globind_flat);
  spcG = spc(globind_flat);
  spvG = spv(globind_flat);
  clear globind_flat

  clear spr spc spv

  fprintf('\nstart dimred ...\n');
  fprintf('N=%d, fill rate %0.0004f (nnz as in a dense nxn matrix with n=%0.0f) ...', ...
          N, ...
          2*(length(sprL)+length(sprG))/N^2, ...
          sqrt(2*(length(sprL)+length(sprG))));


  X = sptsne_optimize(N, L, G, ...
                      sprL, spcL, spvL, ... 
                      sprG, spcG, spvG, ... 
                      tsneopt);
  