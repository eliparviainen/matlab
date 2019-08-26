function sigma2 = apprsigmas(nnd2, u, tol)
%
%  sigma2 = apprsigmas(nnd2, u, tol)
% 
%     Determines kernel widths for t-SNE neighborhood probabilities,
%     using a limited set of pairwise distances (those to nearest
%     neighbors). Binary search for a sigma that makes log(u) equal
%     to the entropy of the kernel, within some tolerance.
%     The code is a straightforward modification of Laurens van der 
%     Maaten's d2p.m code.
% 
%     nnd2           NxL, matrix of squared euclidean distances to L 
%                    nearest neighbors. L should be larger than u.
%     u              neighborhood size as perplexity 
%                    (aka effective number of neighbors)
%     tol (optional) search tolerance
% 
%     output:        kernel width (sigma squared)
%
%
%
% (c) Eli Parviainen, 2014
% Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
% 

    
    if ~exist('tol', 'var') || isempty(tol)
        tol = 1e-4; 
    end
    
    
    % Initialize some variables
    N = size(nnd2, 1);    % number of instances
    beta = ones(N, 1);    % empty precision vector
    logU = log(u);        % log of perplexity (= entropy)

    

    % Run over all datapoints
    for row=1:N
      
        if mod(row, 50)==0
          fprintf('\rrow=%d',row);
        end
        
        % Set minimum and maximum values for precision
        betamin = -Inf; 
        betamax = Inf;

        % Compute the Gaussian kernel and entropy for the current precision

        [H, thisP] = Hbeta(nnd2(row,:), beta(row));
        
        %nonz
        %error
        
        % Evaluate whether the perplexity is within tolerance
        Hdiff = H - logU;
        tries = 0;
        while abs(Hdiff) > tol && tries < 50
            
            % If not, increase or decrease precision
            if Hdiff > 0
                betamin = beta(row);
                if isinf(betamax)
                    beta(row) = beta(row) * 2;
                else
                    beta(row) = (beta(row) + betamax) / 2;
                end
            else
                betamax = beta(row);
                if isinf(betamin) 
                    beta(row) = beta(row) / 2;
                else
                    beta(row) = (beta(row) + betamin) / 2;
                end
            end
            
            % Recompute the values
            [H, thisP] = Hbeta(nnd2(row,:), beta(row));
            Hdiff = H - logU;
            tries = tries + 1;
        end
       
    end    
     
    disp(['Mean value of sigma: ' num2str(mean(sqrt(1 ./ beta)))]);
    disp(['Minimum value of sigma: ' num2str(min(sqrt(1 ./ beta)))]);
    disp(['Maximum value of sigma: ' num2str(max(sqrt(1 ./ beta)))]);

sigma2 = 1 ./ beta;
    
end
    



% Function that computes the Gaussian kernel values given a vector of
% squared Euclidean distances, and the precision of the Gaussian kernel.
% The function also computes the perplexity of the distribution.
function [H, P] = Hbeta(nnd2, beta)

    P = exp(-nnd2 * beta);
    sumP = sum(P);
    H = log(sumP) + beta * sum(nnd2 .* P ) / sumP;
    P = P / sumP;
    
end

