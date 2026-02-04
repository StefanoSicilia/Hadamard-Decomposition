function [W1,H1,W2,H2,info]=HadDec(X,r,opts)
%% HadDec: rank-r Hadamard (entrywise) decomposition of a matrix
% Given an m-by-n matrix X and a rank r, it computes an approximate 
% Hadamard Decomposition X~= (W1*H1').*(W2*H2'), where W1 and W2 are m-by-r
% matrices and H1 and H2 are n-by-r matrices.
% It minimizes the objective function 
% 0.5*norm(X-(W1*H1').*(W2*H2'))^2/norm(X,'fro')^2 
% in three possible ways, depending on the method chosen (see below).
%
% Inputs:
%   X: m-by-n matrix to be decomposed
%   r: a postive integer for the rank-r Hadamard Decomposition
%   opts: struct with fields
%       1) method - the method chosen (see below) [default Manopt]
%       2) init - initialization method (see below) [default all]
%       3) maxit - maximum number of iterations [default 1e6]
%       4) maxtime - maximum time allowed for the computation (this does 
%       not include the error computation for opts.method='manBCDsparse')
%       [default 10 seconds]
%       5) tol - tolerance on the objective function [default 0.5e-30]
%       6) momentum - parameters for the extrapolation 
%       [default [0.75,1,1.05,1.01,1.5]]
%       7) tau (only for manBCDs) - parameter for gradient descent stepsize
%       [default 1.5]
%       8) Wblock (only for manBCDs) - number of iteration per each W block
%       [default 1]
%       9) Hblock (only for manBCDs) - number of iteration per each H block
%       [default 1]
% See the tests for some practical choices of the parameters. Otherwise the
% function can be called by just specifying the rank r as HadDec(X,r) with 
% default parameters.
%
% Outputs:
%   the matrices W1,H1,W2 and H2 of the decomposition
%   a struct 'info' with fields
%       1) err - vector with the relative errors (objective function)
%       2) times - vector with the time required by each iteration
%       3) init - string about the initialization method performed
%
% The methods available are
%   1) manBCD: it uses a 2 block coordinate descend algorithm for the
%   rank-r^2 matrices W=face_split(W1,W2) and H=face_split(H1,H2) and it 
%   optimizes on Bmr and Bnr respectively, where Bmr is the manifold of 
%   matrices that admit a face-split decomposition.
%   2) BCD:it uses the 4 block coordinate descent algorithm described in
%   the paper by Wertz et al., with extrapolation (momentum) provided in 
%   opts.momentum.
%   3) Manopt: it uses a product manifold for X1=W1*H1' and X2=W2*H2' and 
%   it optimizes through Manopt to find an approximation X~=X1.*X2.
%   4) manBCDsparse: like manBCD, but for large sparse matrices.
%
% The initializations available are:
%   1) 'SVD-based' (see Init_SVDbased)
%   2) 'FS' (see Init_FS)
%   3) 'FSL' (see Init_FSL)
%   4) 'FSR' (see Init_FSR)
%   5) 'best' performs the four previous initializations and runs the
%   algorithm only for the one with lowest initial error (see HadDec_init)
%   6) 'all' runs the algorithm for all the first four initializations and
%   selects the best result in terms of final error (see HadDec_all_init)
%   7) 'given' uses as starting points given matrices W1,H1,W2 and H2 
%   8) 'rand' random initialization

    % Default parameters
    if nargin<2
        error('Not enough input arguments: missing matrix and/or desired rank.')
    end

    if nargin>3
        error('Too many input arguments.')
    end

    if nargin==2
        opts=struct('method','Manopt','init','all','maxit',1e6,...
            'maxtime',10,'tol',0.5*1e-30,'tau',1.5,'Hblock',1,'Wblock',1);
        opts.momentum=[0.75,1,1.05,1.01,1.5];
    else 
        % nargin==3
        if ~isfield(opts,'method')
            opts.method='Manopt';
        end
        if ~isfield(opts,'init')
            opts.init='all';
        end
        if ~isfield(opts,'maxit')
            opts.maxit=1e6;
        end
        if ~isfield(opts,'maxtime')
            opts.maxtime=10;
        end
        if ~isfield(opts,'tol')
            opts.tol=0.5*1e-30;
        end
        if ~isfield(opts,'tau')
            opts.tau=1.5;
        end
        if ~isfield(opts,'Wblock')
            opts.Wblock=1;
        end
        if ~isfield(opts,'Hblock')
            opts.Hblock=1;
        end
        if ~isfield(opts,'momentum')
            opts.momentum=[0.75,1,1.05,1.01,1.5];
        end
    end

    % Normalization of X
    normX=norm(X,'fro');
    X=X/normX;

    % Initialization and parameters
    info_init=opts.init;
    switch info_init
        case 'all'
            allinit=HadDec_all_init(X,r);
            opts.maxtime=opts.maxtime/4;
        case 'best'
            [W1,H1,W2,H2,info_init]=HadDec_init(X,r);
        case 'SVD-based'
            [W1,H1,W2,H2]=Init_SVDbased(X,r);
        case 'FS'
            [W1,H1,W2,H2]=Init_FS(X,r);
        case 'FSL'
            [W1,H1,W2,H2]=Init_FSL(X,r);
        case 'FSR'
            [W1,H1,W2,H2]=Init_FSR(X,r);
        case 'given'
            W1=opts.W1;
            H1=opts.H1;
            W2=opts.W2;
            H2=opts.H2;   
            W1=W1/sqrt(normX);
            H2=H2/sqrt(normX);
        case 'rand'
            [m,n]=size(X);
            rng(6)
            W1=rand(m,r);
            H1=rand(n,r);
            W2=rand(m,r);
            H2=rand(n,r);
        case 'sparse'
            allinit=HadDec_init_sparse(X,r);
            opts.maxtime=opts.maxtime/4;
        otherwise
            error('Initialization method not available.')
    end


    % Selection of the method: manBCD, BCD or Manopt
    switch opts.method
        case 'manBCD'
            loop=@(W1,H1,W2,H2) loop_manBCD(X,W1,H1,W2,H2,opts);
        case 'BCD'
            opts.tol=sqrt(2*opts.tol);
            loop=@(W1,H1,W2,H2) loop_BCD(X,W1,H1,W2,H2,opts);
        case 'Manopt'
            loop=@(W1,H1,W2,H2) loop_Manopt(X,W1,H1,W2,H2,opts);
        case 'manBCDsparse'
            loop=@(W1,H1,W2,H2) loop_manBCD_sparse(X,W1,H1,W2,H2,opts);
        otherwise
            error('Method not available')
    end

    

    % Main computation
    if ~strcmp(opts.init,'all') && ~strcmp(opts.init,'sparse')
        % Single inizialization 
        [W1,H1,W2,H2,info]=loop(W1,H1,W2,H2);
        info.init=info_init;
    else
        % Multiple initialization
        W1=allinit.W1_svd;
        W2=allinit.W2_svd;
        H1=allinit.H1_svd;
        H2=allinit.H2_svd;
        [W1_svd,H1_svd,W2_svd,H2_svd,info_svd]=loop(W1,H1,W2,H2);
        W1=allinit.W1_FS;
        W2=allinit.W2_FS;
        H1=allinit.H1_FS;
        H2=allinit.H2_FS;
        [W1_FS,H1_FS,W2_FS,H2_FS,info_FS]=loop(W1,H1,W2,H2);
        W1=allinit.W1_FSL;
        W2=allinit.W2_FSL;
        H1=allinit.H1_FSL;
        H2=allinit.H2_FSL;
        [W1_FSL,H1_FSL,W2_FSL,H2_FSL,info_FSL]=loop(W1,H1,W2,H2);
        W1=allinit.W1_FSR;
        W2=allinit.W2_FSR;
        H1=allinit.H1_FSR;
        H2=allinit.H2_FSR;
        [W1_FSR,H1_FSR,W2_FSR,H2_FSR,info_FSR]=loop(W1,H1,W2,H2);
        err_vec=[info_svd.err(end),info_FS.err(end),...
            info_FSL.err(end),info_FSR.err(end)];
        method={'svd','FS','FSL','FSR'};
        W1_vec={W1_svd,W1_FS,W1_FSL,W1_FSR};
        W2_vec={W2_svd,W2_FS,W2_FSL,W2_FSR};
        H1_vec={H1_svd,H1_FS,H1_FSL,H1_FSR};
        H2_vec={H2_svd,H2_FS,H2_FSL,H2_FSR};
        info_vec={info_svd,info_FS,info_FSL,info_FSR};
        [~,index]=min(err_vec);
        info=info_vec{index};
        info.init=method{index};
        W1=W1_vec{index};
        W2=W2_vec{index};
        H1=H1_vec{index};
        H2=H2_vec{index};
    end
    
    % Output
    alpha=sqrt(sqrt(normX));
    W1=alpha*W1;
    H1=alpha*H1;
    W2=alpha*W2;
    H2=alpha*H2;

end



