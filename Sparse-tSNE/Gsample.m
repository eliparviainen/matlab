function [rows cols] = Gsample(N, G, verbose, onlyUpperTri)

% [rows cols] = Gsample(N, G, verbose, onlyUpperTri)
%
%     Samples random global links for sparse t-SNE. 
%
%     N                        matrix size is NxN
%     G                        number of global links per row
%     verbose (optional)       0=silent (default), 1=report progress
%     onlyUpperTri (optional)  0=both triangles, 1=upper triangle (default)
%                              (full matrix is generated as an intermediate 
%                              result, so 0 saves some effort if you need 
%                              the full matrix later)
%     Results are in format suitable for creating sparse matrices:
%     sp = sparse(rows, cols, 1, N, N);
%
%
% (c) Eli Parviainen, 2014
% Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
% 


if ~exist('onlyUpperTri','var')
  onlyUpperTri = 1;
end;


if ~exist('verbose','var')
  verbose=0;
end;

if verbose, fprintf('\nentering Gsample'); end;

% make even
G=2*ceil(G/2);

rows = uint32(zeros(N,G));
cols = uint32(zeros(N,G));

if verbose, fprintf('\n    memory allocated'); end;

% circulant, diag=0, G/2-bands above+below diag
% (any other well-form starting matrix would do also)
sind=0;
for rind=1:N    
  rows(rind, 1:G)=rind;  
  cols(rind, 1:G/2)=mod(rind-G/2-1:rind-2,N)+1;
  cols(rind, G/2+1:G)=mod(rind:rind+G/2-1,N)+1;
end;

if verbose, fprintf('\n    initial matrix created'); end;

% NxG to 1xN*Gx1
rows=rows(:)';
cols=cols(:)';

% random symmetry-preserving permutation
[~, jarj] = sort(rand(1,N));
rows=jarj(rows);
cols=jarj(cols);

if verbose, fprintf('\n    permutation ready'); end;

if onlyUpperTri
  
  upperTri = find(cols>rows);
  rows = rows(upperTri);
  cols = cols(upperTri);
  
  if verbose, fprintf('\n    upper tri ready'); end;
end; % if upper tri

rows = rows(:);
cols = cols(:);


if verbose, fprintf('\nGsample exiting.'); end;