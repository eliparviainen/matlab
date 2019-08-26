function X = sptsne_optimize(N, numL, numG, ...
                             ProwsUL, PcolsUL, PvalsUL, ...
                             ProwsUG, PcolsUG, PvalsUG, ...
                             opt, oldX)
  %  
  %
  %  X = sptsne_optimize(N, numL, numG, ...
  %                      ProwsUL, PcolsUL, PvalsUL, ...
  %                      ProwsUG, PcolsUG, PvalsUG, ...
  %                      opt, oldX)
  % 
  %   Computes an embedding with sparse t-SNE. Called from sptsne.m.
  %   Uses code snippets from Laurens van der Maaten's tsne_p.m
  %   (notably, same optimizer parameters are used).
  % 
  %   N                number of data points
  %   numL, numG       numbers of local and global links, per row
  %   ProwsUL, PcolsUL, PvalsUL
  %                    local links of the sparse similarity matrix
  %   ProwsUG, PcolsUG, PvalsUG
  %                    global links of the sparse similarity matrix
  %   opt              optimizer settings, see sptsne.m
  %   oldX (optional)  initial configuration of the embedding
  % 
  %
  % (c) Eli Parviainen, 2014
  % Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
  % 
  
  
  % ====================
  % init 
  % ====================
  
  % normalize whole matrix sum to one
  nf = 2*sum(PvalsUL(:))+2*sum(PvalsUG(:));
  PvalsUL = PvalsUL ./ nf;
  PvalsUG = PvalsUG ./ nf;


  ProwsUG=ProwsUG(:); PcolsUG=PcolsUG(:); PvalsUG=PvalsUG(:);
  ProwsUL=ProwsUL(:); PcolsUL=PcolsUL(:); PvalsUL=PvalsUL(:);

  if opt.useEarlyExaggeration
    PvalsUL = PvalsUL * opt.exaggFactor;
    PvalsUG = PvalsUG * opt.exaggFactor;
  end;

  if exist('oldX','var') & ~isempty(oldX)
    X=oldX;
  else
    X = 0.0001 * randn(N, opt.outputdims);
  end;

  y_incs  = zeros(size(X));
  gains = ones(size(X));

  n=size(X,1);

  t0=cputime; 

  % ====================
  % main loop
  % ====================

  for iter=1:opt.max_iter
    
    % ------------------------------
    % current neigborhood probs in low-dim space
    % ------------------------------

    % squared euclidean dists, using (a-b)^2 = a^2 + 2ab + b^2
    sum_X = sum(X .^ 2, 2);      
    sqdistL = -2*sum(X(ProwsUL, :) .* X(PcolsUL, :),2);
    sqdistL = sqdistL + sum_X(PcolsUL) + sum_X(ProwsUL);
    sqdistG = -2*sum(X(ProwsUG, :) .* X(PcolsUG, :),2);
    sqdistG = sqdistG + sum_X(PcolsUG) + sum_X(ProwsUG);
    clear sum_X;
    
    % t-distribution
    tdistUnnormL = 1 ./ (1+sqdistL);
    tdistUnnormG = 1 ./ (1+sqdistG);
    
    % lengths of G-links
    e = sqrt(sqdistG);

    
    % ------------------------------
    % how many points a bubble represents
    % ------------------------------
    
    % -1=self
    bubblerepr = (N-numL-1)/numG;
    
    % estimate current data cloud area in 2D
    Rwhole = sqrt(max(max(sqdistL),max(sqdistG)))/2;
    if isempty(Rwhole) % only occurs with empty L
      Rwhole = sqrt(max(sqdistG))/2;
    end;
    
    % rowG = num bubbles
    R = Rwhole/nthroot(numG,opt.outputdims);
    
    % -1 since bubble center is special case
    wei = bubblerepr-1;
    
    clear sqdistL sqdistG

    % ------------------------------
    % exact gradient for local links
    % ------------------------------
    
    LvalsL_PT = PvalsUL.*tdistUnnormL; % attraction
    LvalsL_TT = -tdistUnnormL.*tdistUnnormL; % repulsion
    
    
    % ------------------------------
    % approximate normalization factor
    % ------------------------------
    
    % 2 since links are upper triangle only
    Z = 2*( sum(tdistUnnormL)+sum( bubblerepr.*tdistUnnormG, 1) );
    % with no G-links would divide by zero without this
    if isempty(Z), Z=1; end;

    clear tdistUnnormL
    
    
    % ------------------------------
    % exact value at bubble center
    % ------------------------------
    
    TTcenter = tdistUnnormG.^2;  clear tdistUnnormG
    
    % ------------------------------
    % average force from rest of the bubble
    % ------------------------------
    
    c1 = e.^2+R.^2;
    c2 = (e-R).^2;
    c3 = (e+R).^2;
    
    a1 = -4.*e.*R./c2;
    a2 = 4.*e.*R./c3;
    
    % always a1<=0, scaled to [0,1]
    A = -a1./(1-a1);
    
    [auxK, auxE] = ellipke(A); 
    K1 = 1./sqrt(1-a1).*auxK;
    clear auxK
    E1 = sqrt(1-a1).*auxE;
    clear auxE
    clear a1 A 

    % if e is close to R, a2>1 can happen (a2 is very close to 1 but larger),
    % prevent it (NaN ok, will be overwritten)
    a2(a2>=1)=NaN;

    [K2, E2] = ellipke(a2); 
    clear a2
    
    I = c2 .* (c1 .* E1 - c3 .* K1) ./ (3 * e .* abs(e-R)) + ...
        c3 .* (c1 .* E2 - c2 .* K2) ./ (3 * e .* (e+R));
    clear c1 c2 c3 E1 E2 K1 K2 e
    
    TTfromBubble = I.*TTcenter./(pi.*R.^2);

    % assumes all Inf's mean e==R
    % (testing for e==R is not enough, since for "e almost R" the 
    % test would fail but K1 would already go to Inf)
    eEqualsR = find(isinf(I)|isnan(I));
    clear I

    TTfromBubble(eEqualsR) = 8/(3*pi)*TTcenter(eEqualsR);
    clear eEqualsR 

    % point self + surrounding bubble
    LvalsG = - wei.*TTfromBubble - TTcenter;
    clear TTcenter TTfromBubble wei

    
    % ------------------------------
    % collect results
    % ------------------------------
    
    LvalsL_QT = LvalsL_TT/Z; clear LvalsL_TT
    LvalsL = LvalsL_PT + LvalsL_QT;  clear LvalsL_PT LvalsL_QT
    LvalsG = LvalsG/Z;
    
    LL = sparse2(ProwsUL, PcolsUL, LvalsL, n, n, length(ProwsUL));  clear LvalsL
    LG = sparse2(ProwsUG, PcolsUG, LvalsG, n, n, length(ProwsUG));  clear LvalsG
    
    % assumes L has zero diagonal
    DL = diag(sum(LL,1)+sum(LL,2)');
    DG = diag(sum(LG,1)+sum(LG,2)');
    
    % ------------------------------
    % (approximate) gradient
    % ------------------------------
    
    
    y_grads = 4 * (DL*X - LL*X - (X'*LL)') + ...
              4 * (DG*X - LG*X - (X'*LG)');
    clear DL DG LL LG
    

    % ------------------------------
    % update the solution
    % ------------------------------    

    % note that the y_grads are actually -y_grads
    gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...         
            + (gains * .8) .* (sign(y_grads) == sign(y_incs));
    gains(gains < opt.min_gain) = opt.min_gain;
    y_incs = opt.momentum * y_incs - opt.epsilon * (gains .* y_grads);
    
    X = X + y_incs;

    % zero mean 
    X = bsxfun(@minus, X, mean(X, 1));

    % ------------------------------
    % optimizer settings
    % ------------------------------    

    
    % Update the momentum if necessary
    if iter == opt.mom_switch_iter
      momentum = opt.final_momentum;
    end
    
    if opt.useEarlyExaggeration && iter == opt.stop_lying_iter
      PvalsUL = PvalsUL ./ opt.exaggFactor;
      PvalsUG = PvalsUG ./ opt.exaggFactor;
    end
    
    if opt.verbose
      % Print out progress
      if ~rem(iter, opt.verboseiter)
        fprintf('\nIteration %d (%s)', iter, str_time(cputime-t0));
        t0=cputime;
      end;
    end
    
  end % main loop
  

% ==============================
% readable time formatting
% ==============================

function str = str_time(t)

  if t<60
    str = sprintf('%0.0f s',t);
  elseif t<3600
    str = sprintf('%0.0f min',t/60);
  else
    str = sprintf('%d h %0.0f min',floor(t/3600),(t-3600*floor(t/3600))/60);
  end;