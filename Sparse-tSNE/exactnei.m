function [nei_ind nei_sqd] = exactnei(x, K, opt, konly, verbose)

% EXACTNEI(x, K, opt, konly)
%    Finds 1st to Kth neighbors for each point in x (NxP, N samples of
%    dimension P), or the Kth neighbor (if konly is set). Opt can be
%    empty or contain field 'block' which determines how many rows are
%    handled at once (saves memory with large x, default is all
%    rows). Returns neighbor indices in nei_ind (size NxK or Nx1)
%    and squared euclidean distances to them in nei_sqd (same
%    size). Verbose can be 0 (no output) or 1 (progress indicated).
%
%    See also: APPRNEI
%
%
% (c) Eli Parviainen, 2014
% Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
% 


if  ~exist('verbose','var')
  verbose = 1;
end;

if  ~exist('konly','var') | isempty(konly)
  konly=0;
end;

if  ~exist('opt','var') | isempty(opt)
  opt=struct();
  % lower memory needs by handling a block of rows
  % instead of all at once 
  opt.block = size(x,1);
end;


N = size(x,1);

sqself = sum(x .^ 2, 2);  


if konly
  % if only one is needed, save memory
  nei_ind = uint32(zeros(N, 1));
  if nargout>1
    nei_sqd = zeros(N, 1);
  end;
  
else
  nei_ind = uint32(zeros(N, K));
  if nargout>1
    nei_sqd = zeros(N, K);
  end;
end;

% loop over data points
if verbose, fprintf('\n'); end;
rangeA=1; 
rangeL=min(opt.block,N);
for ind=1:ceil(N/opt.block)

  if mod(ind,10)==0 & verbose
    fprintf('b');
  end;
  
  % squared distance from block points to all points
  sqd=bsxfun(@plus,bsxfun(@plus,-2*x(rangeA:rangeL,:)*x',...
                          sqself(rangeA:rangeL)),sqself');

  % sort is col-wise, make it row-wise, sort, then return to col-wise
  [sqd si] = sort(sqd', 'ascend');
  
  % pick K nearest
  si = si(1:K, :)';
  sqd = sqd(1:K, :)';

  if konly
    si = si(:, K);
    sqd = sqd(:, K);
  end;


  nei_ind(rangeA:rangeL, :) = uint32(si);
  clear si
  if nargout>1
    nei_sqd(rangeA:rangeL, :) = sqd;
    clear sqd
  end;
  
  rangeA=rangeA+opt.block;
  rangeL=rangeL+opt.block;
  rangeL=min(rangeL, N);


end;
if verbose, fprintf('\n'); end;
