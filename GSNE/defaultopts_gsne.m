function gsneopt = defaultopts_gsne(numG, walk_stop_thr, local_nei_thr, numiter, verbose)
%
% gsneopt = defaultopts_gsne(numG, walk_stop_thr, local_nei_thr, maxiter, verbose)
%
%    Fills the critical entries of gsne options with user-supplied values, 
%    and other fields with default values.
%
%    numG                Number of global links. First try 5..50, increase if needed.
%    walk_stop_thr       Neighborhood size parameter (typical values 1/10, 1/20, 1/50).
%    local_nei_thr       Size of locally linked area (e.g. 0.01 to 0.1 times walk_stop_thr).
%    numiter (optional)  Number of iterations, default 500.
%    verbose (optional)  0=silent, 1=report progress (default).
%
%    Defaulted fields include 
%          walk_maxsteps = 100   (good expanders might need smaller value?)  
%          useweights = 0        (set to 1 for weighted graphs) 
%          lowmem = 0            (set to 1 for slower but memory-efficient computation 
%                                 of neighborhoods)
%
% (c) Eli Parviainen, 2014
% Use FREELY for any NON-COMMERCIAL purpose, at your OWN RISK.
% 

  if ~exist('verbose','var')
    verbose = 1;
  end;


  if ~exist('numiter','var')
    numiter = 500;
  end;

  
  % Most of these are internal optimizer parameters of 
  % sparse t-SNE, which seldom need tuning.
  tsneopt = struct(); 
  tsneopt.max_iter = numiter;
  tsneopt.outputdims = 2;
  tsneopt.momentum = 0.5000;
  tsneopt.final_momentum = 0.8000;
  tsneopt.mom_switch_iter = 250;
  tsneopt.useEarlyExaggeration = 1;
  tsneopt.stop_lying_iter = 100;
  tsneopt.exaggFactor = 4;
  tsneopt.epsilon = 500;
  tsneopt.min_gain = 0.0100;
  tsneopt.plotprogress = 0;
  tsneopt.verbose = verbose;
  tsneopt.verboseiter = 100;


  % GSNE options
  gsneopt = struct(); 
  
  gsneopt.tsneopt = tsneopt;
  gsneopt.verbose = verbose;
  
  % matrix of probabilities is determined with
  % successive matrix powers if lowmem=0, and with
  % a slower but memory-saving for-loop if lowmem=1
  gsneopt.lowmem = 0; 

  % with 1, use whatever is in data without any checking
  % (read: YOU must make sure edge weights are similarities)
  % with 0, all edges are made to have weight 1
  gsneopt.useweights = 0;

  % number of global links
  gsneopt.numG = numG;

  % stop the lazy random walk when the starting node
  % has this little probability mass left
  gsneopt.walk_stop_thr = walk_stop_thr;
  
  % safeguard the walk from becoming overly long
  gsneopt.walk_maxsteps = 100;

  % any node reachable with higher probability is 
  % a local neighbor (L-linked node)
  gsneopt.local_nei_thr = local_nei_thr;

  