function X = gsne_oldsparse(N, Lspr, Lspc, Lspv, gsneopt, verbose)
  
%  -----------------------------------------------------------
%  This is a version of gsne.m that uses original 
%  Matlab sparse. The version with SuitSparse sparse2 function
%  is faster and needs less memory.
%  -----------------------------------------------------------
%
% X = gsne(N, Lspr, Lspc, Lspv, gsneopt, verbose)
%
%    Computes a graph embedding with GSNE. 
%
%    N                  number of nodes
%    Lspr, Lspc, Lspv   graph in sparse matrix format, upper triangle only
%    gsneopt            options, see defaultopt_gsne.m
%    verbose (optional) 0=silent, 1=report progress
%    
%    X                  Nx2 matrix, node coordinates 
%                       usage example: plot(X(:,1), X(:,2), '.')
%
%
% (c) Eli Parviainen, 2014
% Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
% 

  
% | off since octave does not short-circuit
  if ~exist('verbose','var') %| isempty(verbose)
    verbose = 1;
  end;

  % not a good idea to use defaults for the critical parameters,
  % but at least the code can be tested
  if ~exist('gsneopt','var') %| isempty(gsneopt)
    warning(sprintf(['\nGSNE: No options given, using arbitrary ' ...
                     'fixed values\n(and these DO affect the results)!']));
    gsneopt = defaultopts_gsne(20, 1/10, 1/1000, 500, verbose);
  end;

  verbose = gsneopt.verbose;  
      
  % clean up duplicates after ...N nodes or ...ITER iterations
  LOWMEM_CLEANN = 1000;
  HIMEM_CLEANITER = 20;
  
  % report progress each ...N nodes
  LOWMEM_VERBOSEN = 1000;


  % =========================================================
  % global links
  % =========================================================


  if verbose, fprintf('\nsample %d global links ... ',gsneopt.numG); end;

  % sample random global links, keep both upper and lower triangle
  upperonly = 0;
  Gverbose = 0; % too much detail for normal use, for debug only
  [Gspr Gspc] = Gsample(N, gsneopt.numG, Gverbose, upperonly);
  clear Gverbose upperonly

  % =========================================================
  % local neighborhoods and probabilities from a random walk
  % =========================================================

  if verbose, fprintf('\nlocal links and neighborhood probs '); end;

  % remove any weights if an unweighted graph was called for
  if ~gsneopt.useweights, Lspv(:) = 1; end;

  
  % ----------------------------------------
  % transition matrix
  % ----------------------------------------
  
  Lspd = (1:N)';

  % Diagonal is filled with rows sums
  rsums = zeros(N,1);
  for ind=1:length(Lspr)
    % assumes Lspr, Lspc specify an upper triangular matrix
    rsums(Lspr(ind)) = rsums(Lspr(ind)) + Lspv(ind);
    rsums(Lspc(ind)) = rsums(Lspc(ind)) + Lspv(ind);
  end;
  clear ind

  % Make rows sum to one.
  % Now, row k has total mass on 2*rsums[k], half at diagonal.
  % After normalizing, diagonal will be 1/2.
  Lspvv = [Lspv; Lspv; rsums] ./ (2*rsums([Lspr; Lspc; Lspd]));
  clear rsums

  % Markov transition matrix
  % MATLAB SPARSE
  sp = sparse(double([Lspr; Lspc; Lspd]), double([Lspc; Lspr; Lspd]), Lspvv, N, N);
  clear Lspr Lspc Lspd Lspvv

  if isfield(gsneopt,'lowmem') & gsneopt.lowmem

    % ----------------------------------------
    % slower, saves memory
    % ----------------------------------------

      if verbose, fprintf('(with the low memory version) ...\n    '); end;
    
    Lspr=[]; Lspc=[]; Lspv=[];
    Gspv = zeros(size(Gspr));
    
    for nind = 1:N        
      
      if mod(nind,LOWMEM_VERBOSEN)==0 & verbose
        fprintf('.');
      end;
      
      % start a walk from node nind
      v = zeros(1, N); v(nind)=1;                     
      ctr=0;
      
      % take steps until mass at node falls below threshold or maximum steps
      while v(nind)>=gsneopt.walk_stop_thr & ...
            ctr<gsneopt.walk_maxsteps
        v=v*sp;
        ctr = ctr+1;
      end;        

      % locations of G-links are known, values are recorded here
      Growind=find(Gspr==nind);  
      Gspv(Growind) = v(Gspc(Growind));
      clear Growind               
            
      % node indices, can index rows and columns
      newind = find(v>gsneopt.local_nei_thr);

      % directed links
      Lspr = [Lspr; ones(length(newind),1)*nind;];
      Lspc = [Lspc; newind(:)];
      Lspv = [Lspv; v(newind)'];
      
      % clean up memory, typically many duplicates
      if mod(nind,LOWMEM_CLEANN)==0 
        if verbose, fprintf('c'); end;

        % same duplicate removal algorithm as elsewhere
        % (a more detailed commentary below)
        Lind = (Lspc-1)*N+Lspr;
        [Lind si] = sort(Lind,'ascend');
        Lspv = Lspv(si);
        notdup = find([0; Lind(1:end-1)]~=Lind);
        Lind = Lind(notdup);
        Lspv = Lspv(notdup);
        clear notdup        
        Lspc = ceil(double(Lind)/double(N));
        Lspr = double(Lind)-(Lspc-1)*double(N);
        clear Lind
      end; % if cleanupiter
      
    end; % for n
    clear v
    clear newind notdup auxv
    clear nind
    
  else
    % ----------------------------------------
    % faster, needs more memory
    % ----------------------------------------

    if verbose, fprintf('(with the matrix power version) ...\n    '); end;
    
    % no local links yet
    Lspr=[]; Lspc=[]; Lspv=[];
    
    % L-links are determined during the walk, but G-links must
    % be read from Gspr, Gspc and values recorded. 
    % Values in <mask> at R,C is location of (R,C) in vector Gspv.    
    % MATLAB SPARSE
    mask = sparse(double(Gspr), double(Gspc), 1:length(Gspr), N, N);
    
    % Short walk may not traverse all global links. Their
    % values become "small".
    Gspv = zeros(length(Gspr),1)+eps;
    

    % start a walk from each node (all probmass at diagonal)
    V = speye(N); 
    
    % names of remaining nodes
    diagC = (1:N)'; 
    
    % how many nodes remain
    M=N; 
    
    % loop until all nodes go below threshold or maximum steps
    ctr=1;
    while M>0 & ctr<gsneopt.walk_maxsteps
      
      if mod(ctr,5)==0 & verbose, fprintf('.'); end;        

      % random walk step
      V=V*sp;   
      
      % Indices of "diagonal" entries of V, i.e. those that correspond to diagonal
      % of the square sp. V is square in the beginning but becomes rectangular 
      % when nodes are dropped. Column location <diagC> of the "diagonal" is in same
      % as before but row location changes every iteration, so the index changes.      
      diagind = (diagC-1)*M+(1:M)'; 

      % Remaining probmass at these nodes ("diagonal" of V monitors this, lazy walk)
      % fell below threshold. The values are saved, and from the next iteration,
      % these nodes are not incluced in V anymore (their random walks stop).
      underthr = find(V(diagind)<gsneopt.walk_stop_thr);
     
      % underC is already a valid index for the matrix sp
      [underR underC] = find(V(underthr,:)>(gsneopt.local_nei_thr));
      
      % underR cannot index V, but only its underthr-subblock
      underR = underthr(underR);      
      % now underR is indices to this-iteration-V, which has the
      % same size as diagC
      
      % rows+cols to index      
      % note that we want to index full V, so there are size(V,1) rows, 
      % not N or length(underthr) 
      underI = (underC-1)*size(V,1)+underR;
      underV = full(V(underI)); clear underI;
      
      % diagC is the location of diagonal of the current row,
      % i.e. the original row number
      underR = diagC(underR);
      % now underR can index full sp matrix or the original nodes
      
      
      % record links as directed 
      Lspr = [Lspr; underR(:)];
      Lspc = [Lspc; underC(:)];          
      Lspv = [Lspv; underV(:)];    
      
      
      % SNind indexes subblocks of V and mask
      SNind = find(mask(underthr,:));
      
      % this is slow, possible to index directly somehow?
      submask = mask(underthr,:);        
      
      % numbers in the mask are indices to Gspv
      NNind = submask(SNind);  clear submask      
      subV = V(underthr, :);
      Gspv(NNind) = subV(SNind);  clear subV
      
      
      % typically many duplicates, clean up every now
      % and then to save memory
      if mod(ctr,HIMEM_CLEANITER)==0
        if verbose, fprintf('c'); end;
        
        % same duplicate removal algorithm as elsewhere
        % (a more detailed commentary below)
        Lind = (Lspc-1)*N+Lspr;
        clear Lspr Lspc
        [Lind si] = sort(Lind,'ascend');
        Lspv = Lspv(si);
        notdup = find([0; Lind(1:end-1)]~=Lind);
        Lind = Lind(notdup);
        Lspv = Lspv(notdup);        
        Lspc = ceil(double(Lind)/double(N));
        Lspr = double(Lind)-(Lspc-1)*double(N);
        clear Lind notdup
      end; % if cleanup iter

      overthr = find(V(diagind)>=gsneopt.walk_stop_thr);
      oldM=M;
      M=length(overthr);
      diagC = diagC(overthr);
      V=V(overthr,:);
      mask=mask(overthr,:); % G-links needs this
      ctr=ctr+1;        
    end;
    clear overthr diagind underthr 
    clear ctr oldM

    
    % if maxsteps stopped the walk, check the remaining V
    if M>0
      % underthr etc are not needed, use all V;
      % otherwise like above
      
      [underR underC] = find(V>gsneopt.local_nei_thr);
      underI = (underC-1)*size(V,1)+underR; 
      underV = V(underI); clear underI;     
      underR = diagC(underR);
      
      Lspr = [Lspr; underR(:)];
      Lspc = [Lspc; underC(:)];      
      Lspv = [Lspv; underV(:)];
      
      SNind = find(mask);     
      NNind = mask(SNind);
      Gspv(NNind) = V(SNind);
      
      clear NNind SNind mask V
    end;    

    clear underC underR diagC
    clear V 
    clear underV M
    
  end; % if lowmem




  
  % =========================================================
  % clean the similarity matrix
  % =========================================================
  % available now: Gspr, Gspc, Gspv; Lspr, Lspc, Lspv
  % G is both upper and lower triangle, symmetric links, asymmetric values
  % L is both upper and lower triangle, asymmetric links, asymmetric values
  
  % combine L- and G-links but remember which is which
  islocallink = uint8([ones(length(Lspr), 1); zeros(length(Gspr), 1)]);
  LGspr = [Lspr; Gspr]; clear Lspr Gspr
  LGspc = [Lspc; Gspc]; clear Lspc Gspc
  LGspv = [Lspv; Gspv]; clear Lspv Gspv

    
  % ----------------------------------------
  % remove duplicates
  % ----------------------------------------

  if verbose, fprintf('\nremoving L/G overlap ...'); end;
  
  % L has no duplicates (at each row, only (r,c) was saved, not (c,r))
  % G has no duplicates due to its construction
  % L and G can overlap with each other, this removes such
  
  % sort indices to matrix entries (makes duplicates adjacent),
  % non-duplicate is where adjacent values are not the same
  LGind = (LGspc-1)*N+LGspr;
  [LGind si] = sort(LGind,'ascend');
  islocallink = islocallink(si);
  LGspv = LGspv(si);
  clear si
  notdup = find([0; LGind(1:end-1)]~=LGind);
  
  % only keep non-duplicates
  LGind = LGind(notdup);
  islocallink = islocallink(notdup);
  LGspv = LGspv(notdup);
  clear notdup
  
  % return from indices to row,col
  LGspc = ceil(double(LGind)/double(N));
  LGspr = double(LGind)-(LGspc-1)*double(N);
  clear LGind
  
  LGspc = uint32(LGspc);
  LGspr = uint32(LGspr);  

  if verbose, fprintf('\nnormalize, symmetrize, uppertri ...'); end;
  
  % ----------------------------------------
  % row-wise normalization
  % ----------------------------------------

  % make matrix rows sum to one
  rsums = zeros(N,1);
  for ind=1:length(LGspr)
    % assumes LGspr,LGspc has both upper and lower triangle
    rsums(LGspr(ind)) = rsums(LGspr(ind)) + LGspv(ind);
  end;
  LGspv = LGspv./rsums(LGspr);
  clear rsums ind
  

 
  % ----------------------------------------
  % symmetrize
  % ----------------------------------------
   
  % islocallink becomes invalid below since indices can reorder,
  % now negative values record G-link locations
  isglobal = find(~islocallink);
  clear islocallink
  LGspv(isglobal) = -LGspv(isglobal);
  clear isglobal
  
  % If a and b are neighbors (L-linked) from both directions,
  % matrix entry receives two values, which sparse2 sums.
  % Due to /2, result is the average value.
  % If there is only one link between a and b, its strength
  % is halved. This is intentional, relationship is seen as
  % weaker if only one participant advocates it.
  %
  % G-links survive this intact since they are symmetric already.
  %
  % MATLAB SPARSE
  sp = sparse(double([LGspr; LGspc]), double([LGspc; LGspr]), [LGspv; LGspv]/2, N, N);
  clear LGspr LGspc LGspv
  [spr spc spv] = find(sp); clear sp;
  spr = uint32(spr);
  spc = uint32(spc);

  
  % ----------------------------------------
  % limit to the upper triangle
  % ----------------------------------------
     
  upperTri = find(spc>spr);
  spr = spr(upperTri);
  spc = spc(upperTri);
  spv = spv(upperTri);
  clear upperTri

  

  % =========================================================
  % embed the graph
  % =========================================================
  
  if verbose, fprintf('\ngraph embedding...'); end;

  % separate local and global links 
  % globality is marked by negative spv values
  
  %lokind = find(islocallink); globind = find(~islocallink); clear islocallink
  lokind = find(spv>=0);
  sprL = spr(lokind); spcL = spc(lokind); spvL = spv(lokind); clear lokind
  globind = find(spv<0);
  sprG = spr(globind); spcG = spc(globind); spvG = -spv(globind); clear globind
  clear spr spc spv

 
  % this is used only for determining bubble radii, avg value will do
  L = ceil(2*length(sprL)/N);
  % MATLAB SPARSE
  X = sptsne_optimize_oldsparse(N, L, gsneopt.numG, ...
                      sprL, spcL, spvL, sprG, spcG, spvG, ... 
                      gsneopt.tsneopt);

