function [spr spc spw neid2] = lggraph(x, numL, numG, neiopt, verbose)
  
% [spr spc spw neid2] = lggraph(x, numL, numG, neiopt, verbose)
% 
%    Local and global links for sparse t-SNE. Local links are to nearest
%    neighbors, from any direction (so the realized number of links per
%    point can differ from numL). Global links are uniformly random. 
%    Overlap of local and global links is removed. By default, neighbors
%    are approximate (= determined in a lower-dimensional space).
% 
%    IN:
%    x                      data, NxP
%    numL                   numL local links per point 
%    numG                   numG global links per point 
%    neiopt (optional)      
%        neiopt.dim         dimension for neighbor search (default 50)
%        neiopt.block       saves memory by running neighbor search
%                           as blocks, not for all data at once (default N)
%        neiopt.exactdist   if present and true, use exact search (=ignore dim)
%    verbose (optional)     0=silent, 1=report progress (default)
%
%    OUT:
%    spr, spc               link locations (usage: sp=sparse(spr,spc,1,N,N);)
%    spw                    link type, 1=local, 2=global; same size as spr
%    neid2                  NxL, squared distances to numL closest neighbors
%
%
% (c) Eli Parviainen, 2014
% Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
% 
  
  
  if ~exist('verbose','var')
    verbose = 1;
  end;
  
  if ~exist('neiopt','var') | isempty(neiopt)
    neiopt=struct();
    neiopt.dim = 50;
    neiopt.block = size(x,1);
  end;
  
  
  N = size(x,1);

  if isfield(neiopt,'exactdist') & neiopt.exactdist
    if verbose, fprintf('\nL-links, exact neighbors ...'); end;
    [nei_ind neid2] = exactnei(x, numL, neiopt);
    
  else
    if verbose, fprintf('\nL-links, approximate neighbors ...');  end;
    [nei_ind neid2] = apprnei(x, numL, neiopt);
    
  end; 
  clear x
  

  sprL = repmat(1:N, 1, size(nei_ind,2));
  spcL = nei_ind;
  
  % neid2 is only needed for computation of sigmas, 
  % not every caller wants it
  if nargout<4
    clear nei_ind neid2
  end;

  
  if verbose,  fprintf('\nsymmetrize L-links...'); end;

  
  spr = double([sprL(:); spcL(:) ]);
  spc = double([spcL(:); sprL(:) ]);
  sprL = spr; clear spr
  spcL = spc; clear spc
  spwL = uint8(ones(size(sprL))); % mark local links as 1
  
  upperTri = find(spcL>sprL);
  sprL = sprL(upperTri);
  spcL = spcL(upperTri);
  spwL = spwL(upperTri);
  clear upperTri
  

  if verbose, fprintf('\nG-linkit...'); end;
  [sprG spcG] = Gsample(N, numG, 1-verbose);

  
  % mark global links as 2
  spwG = uint8(2*ones(size(sprG)));

  spr = [sprL(:); sprG(:)]; clear sprL sprG
  spc = [spcL(:); spcG(:)]; clear spcL spcG
  spw = [spwL(:); spwG(:)]; clear spwL spwG
  
  
  if verbose, fprintf('\nremove duplicates ...'); end;
  
  % row/col to ind
  LGind = (spc-1)*N+spr;
  clear spr spc
  [LGind,si] = sort(LGind,'ascend');
  spw = spw(si);
  % duplicates are now sequences of same val, choose
  % values with a different neighbor
  keep = find([0; LGind(1:end-1)]~=LGind);
  LGind = LGind(keep);
  spw = spw(keep);
  
  % ind to row/col
  spc = ceil(double(LGind)/double(N));
  spr = double(LGind)-(spc-1)*double(N);
  clear LGind

  
  if verbose, fprintf('\nlggraph exiting.'); end;