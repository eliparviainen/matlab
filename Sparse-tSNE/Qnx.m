function Q = Qnx(x, X, KK, S)
%
% Q = QNX(x, X, KK, S)
%     Value of Qnx quality criterion (average number of shared points
%     in high- and low-dimensional neighborhoods of same size).
%
%     input:
%     x   original data, NxP1
%     X   low-dimensional data, NxP2
%     KK  vector neighborhood sizes, integers between 1 and N-1
%     S   (optional) evaluate Qnx on a subsample of S points, not all points
%
%     output:
% 
% Reference: 
%     John A. Lee and Michel Verleysen, 2009.
%     Quality assessment of dimensionality reduction: Rank-based criteria.
%     Neurocomputing, vol 72, pp. 1431-1443.
%
%
% (c) Eli Parviainen, 2014
% Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
% 


N=size(X,1);

if ~exist('S','var')
  S=N;
end;

% evaluate quality for subsample of given size
subsample = randperm(N);
subsample = subsample(1:min(S,N));

if ~exist('verbose','var')
  verbose=1;
end;

% --------------------------------------------------

% input data neighbors

% dt2 is SxN, % squared distances from x(subsample, :) to x
sq2 = sum(x .^ 2, 2); 
sq1 = sq2(subsample);
dt2=bsxfun(@plus,bsxfun(@plus,-2*x(subsample,:)*x',sq1),sq2');
clear sq1 sq2

% ordered indices to all N neighbors of the S subsample points
[~, neiHi] = sort(dt2','ascend'); clear dt2
neiHi=neiHi';
neiHi = neiHi(:, 2:end); % first point is self, exclude


% --------------------------------------------------

% low dim neighbors 

% dt2 is SxN, % squared distances from X(subsample, :) to X
sq2 = sum(X .^ 2, 2); 
sq1 = sq2(subsample);
dt2=bsxfun(@plus,bsxfun(@plus,-2*X(subsample,:)*X',sq1),sq2');
clear sq1 sq2

% ordered indices to all N neighbors of the S subsample points
[~, neiLo] = sort(dt2','ascend'); clear dt2
neiLo=neiLo';
neiLo = neiLo(:, 2:end); % exclude self


% --------------------------------------------------
if verbose, fprintf('\nQnx: '); end;


% precision of KK-sized neighborhoods of S
prec = zeros(size(neiHi,1), length(KK));

% each scale (neighborhood size) is checked individually
for kind = 1:length(KK)
  K = KK(kind);

  if mod(kind, 50)==0 & verbose
    fprintf('.');
  end;
  
  % loop over points
  for pind = 1:size(neiHi,1)

    % auxind implements a faster version of intersect(relevant, retrieved)
    auxind = zeros(1,N);
    relevant = neiHi(pind, 1:K);
    auxind(relevant)=1; % right point
    
    retrieved = neiLo(pind,1:K);
    auxind(retrieved)=auxind(retrieved)+2; % in right place
    
    % auxind==1+2==3, right point in right place
    relretv = find(auxind==3);
    
    % percentage of neighbors
    prec(pind, kind) = length(relretv)/K;
  end; % for point
end; % for neighborhood size

% --------------------------------------------------
% average to get Qnx-criterion

Q = mean(prec, 1); % mean over points
if verbose, fprintf('\n'); end;