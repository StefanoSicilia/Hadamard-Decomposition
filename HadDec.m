function [W1,H1,W2,H2,info]=HadDec(X,r,opts)
%% HadDec: rank-r Hadamard (entrywise) decomposition of a matrix
% Given an m-by-n matrix X and r=[r1,r2], it computes an approximate 
% Hadamard Decomposition X~=(W1*H1').*(W2*H2'), where W1 and W2 are m-by-r1
% matrices and H1 and H2 are n-by-r2 matrices.
% It minimizes the objective function 
% 0.5*norm(X-(W1*H1').*(W2*H2'))^2/norm(X,'fro')^2 
% in four possible ways, depending on the method chosen (see below).
% See the tests for some practical choices of the parameters. 
%
% The function can be called by just specifying the ranks r=[r1,r2] as 
% [W1,H1,W2,H2,info]=HadDec(X,r) with default parameters (see below), or
% even easier [W1,H1,W2,H2,info]=HadDec(X,r) with r an integer greater than
% 1 such that r1=r2=r.
%
% Inputs:
%   X: m-by-n matrix to be decomposed
%   r: an integer greater than 1 for the rank-r Hadamard Decomposition
%       or a vector with two ranks [r1,r2] greater than 1
%   opts: struct with fields
%       1) method - the method chosen (see below) [default 'Manopt']
%       2) init - initialization method (see below) [default 'FS']
%       3) maxit - maximum number of iterations [default 1e6]
%       4) maxtime - maximum time allowed for the computation (this does 
%       not include the error computation, since it might be removed)
%       [default 10 seconds]
%       5) tol - tolerance on the objective function [default 1e-16]
%       6) momentum (not for Manopt) - parameters for the extrapolation in 
%       the form [beta, beta_tilde, gamma, gamma_tilde, eta]: beta is the 
%       extrapolation parameter, while the rest of the values are used to
%       tune beta during the iteration (upper bounds, decreasing and
%       increasing factors) [default [0.25,1,1.05,1.01,1.5] for manBCD and
%       projBCD, [0.75,1,1.05,1.01,1.5] for BCD]
%       7) tau (only for manBCD and projBCD) - parameter for gradient 
%       descent stepsize [default 0.95]
%       8) Iter_W (only for manBCD and projBCD) - number of iterations per 
%       each W block [default 2]
%       9) Iter_H (only for manBCD and projBCD) - number of iterations per 
%       each H block [default 2]
%       10) sparsity - flag to compute the error differently in the matrix   
%       sparse case [default 0]
%       11) noloops - flag for avoiding implementation loops [default 1]
%       12) rescale - flag for activating the rescaling procedure (see
%       normalize function) [default 1]
%       13) theta - parameter for the perturbation of the zero rows in the
%       initial matrices [default 1e-4]
%
% Outputs:
%   the matrices W1,H1,W2 and H2 of the decomposition
%   a struct 'info' with fields
%       1) err - vector with the relative errors 
%       norm(X-(W1*H1').*(W2*H2'))/norm(X,'fro')
%       2) times - vector with the time required by each iteration
%       3) init - string about the initialization method performed
% and supplementary fields in the case of multiple initializations
%       4) global - info of all the intializations performed
%       5) ratioinit - relative ratio between the best initialization and the
%       other ones performed, namely |e_best-e_init|/(|e_best|+1e-8)
%
% The methods available are [default 'Manopt']:
%   1) Manopt: it uses a product manifold for X1=W1*H1' and X2=W2*H2' and 
%   it optimizes through Manopt to find an approximation X~=X1.*X2.
%   2) manBCD: it uses a 2 block coordinate descend algorithm for the
%   rank-(r1*r2) matrices W=face_split(W1,W2) and H=face_split(H1,H2) and  
%   it optimizes on Bmr and Bnr respectively, where Bmr is the manifold of 
%   matrices that admit a face-split decomposition with ranks r1 and r2.
%   3) projBCD: alternative projections onto the manifold Bmr and Bnr, using
%   the representation X~=WH', with W in Bmr and H in Bnr.
%   4) BCD:it uses the 4 block coordinate descent algorithm described in
%   the paper by Wertz et al., with extrapolation (momentum) provided in 
%   opts.momentum.
%
% The initializations available are [default 'all'] :
%   1) 'SVD-based' (see Init_SVDbased)
%   2) 'FS' (see Init_FS)
%   3) 'FSL' (see Init_FSL)
%   4) 'FSR' (see Init_FSR)
%   5) 'best' performs the four previous initializations and runs the
%   algorithm only for the one with lowest initial error (see HadDec_init)
%   6) 'all' runs the algorithm for all the first four initializations and
%   selects the best result in terms of final error (see HadDec_all_init)
%   7) 'given' uses as starting points given matrices W1=opts.W1, 
%   H1=opts.H1, W2=opts.W2 and H2=opts.H2. 
%   8) 'rand' random initialization

    % Initial checks and default parameters
    if nargin<2
        error('Not enough input arguments: missing matrix and/or desired rank(s).')
    end
    switch length(r)
        case 1
            r1=r;
            r2=r;
        case 2
            r1=r(1);
            r2=r(2);
        otherwise
            error('There must be 1 rank (r) or 2 ranks (r1 and r2).')
    end
    if min([r1,r2])<=1
        error('Each rank must be greater than 1.')
    end
    if nargin>3
        error('Too many input arguments.')
    end
    if nargin==2
        opts=[]; 
    end
    if ~isfield(opts,'method')
        opts.method='Manopt';
    end
    if ~isfield(opts,'init')
        if all(isfield(opts,{'W1','H1','W2','H2'}))
            opts.init='given';
        else
            opts.init='FS';
        end
    end
    if ~isfield(opts,'maxit')
        opts.maxit=1e6;
    end
    if ~isfield(opts,'maxtime')
        opts.maxtime=10;
    end
    if ~isfield(opts,'tol')
        opts.tol=1e-15;
    end
    if ~isfield(opts,'tau')
        opts.tau=0.95;
    end
    if ~isfield(opts,'theta')
        opts.theta=1e-4;
    end
    if ~isfield(opts,'Iter_W')
        opts.Iter_W=2;
    end
    if ~isfield(opts,'Iter_H')
        opts.Iter_H=2;
    end
    if ~isfield(opts,'sparsity')
        opts.sparsity=issparse(X);
    end
    if ~isfield(opts,'noloops')
        opts.noloops=1;
    end
    if ~isfield(opts,'momentum')
        if strcmp(opts.method,'BCD')
            opts.momentum=[0.75,1,1.05,1.01,1.5];
        else
            opts.momentum=[0.25,1,1.05,1.01,1.5];
        end
    end
    if ~isfield(opts,'rescale')
        opts.rescale=1;
    end

    % Normalization of X and zero matrix case
    normX=norm(X,'fro');
    if normX==0
        [m,n]=size(X);
        W1=zeros(m,r1);
        H1=zeros(n,r1);
        W2=zeros(m,r2);
        H2=zeros(n,r2);
        info=struct('err',0,'time',0,'init',[],'global',[],'ratioinit',[]);
        return
    end
    X=X/normX;

    % Initialization and parameters
    info_init=opts.init;
    opt_loop=opts.noloops;
    opt_sparse=opts.sparsity;
    ranks=[r1 r2];
    switch info_init
        case 'all'
            if prod(ranks)>min(size(X))
                error('"all" init. can only be used when r1*r2<=min(m,n)')
            end
            WH=Init_all(X,ranks,opt_loop,opt_sparse);
            opts.maxtime=opts.maxtime/4;
        case 'best'
            [W1,H1,W2,H2,info_init]=Init_best(X,ranks,opt_loop,opt_sparse);
        case 'SVD-based'
            [W1,H1,W2,H2]=Init_SVDbased(X,ranks,opt_sparse);
        case 'FS'
            if prod(ranks)>min(size(X))
                error('FS init. can only be used when r1*r2<=min(m,n)')
            end
            [W1,H1,W2,H2]=Init_FS(X,ranks,opt_loop,opt_sparse);
        case 'FSL'
            if prod(ranks)>min(size(X))
                error('FSL init. can only be used when r1*r2<=min(m,n)')
            end
            [W1,H1,W2,H2]=Init_FSL(X,ranks,opt_loop,opt_sparse);
        case 'FSR'
            if prod(ranks)>min(size(X))
                error('FSR init. can only be used when r1*r2<=min(m,n)')
            end
            [W1,H1,W2,H2]=Init_FSR(X,ranks,opt_loop,opt_sparse);
        case 'given'
            W1=opts.W1;
            H1=opts.H1;
            W2=opts.W2;
            H2=opts.H2;   
            sqsqX=sqrt(sqrt(normX));
            W1=W1/sqsqX;
            H1=H1/sqsqX;
            W2=W2/sqsqX;
            H2=H2/sqsqX;
        case 'rand'
            rng(1)
            [m,n]=size(X);
            W1=rand(m,r1);
            H1=rand(n,r1);
            W2=rand(m,r2);
            H2=rand(n,r2);
        otherwise
            error('Initialization method not available.')
    end

    % Selection of the method: BCD, Manopt, manBCD or proj
    switch opts.method
        case 'Manopt'
            loop=@(W1,H1,W2,H2) loop_Manopt(X,W1,H1,W2,H2,opts);
        case 'manBCD'
            loop=@(W1,H1,W2,H2) loop_manBCD(X,W1,H1,W2,H2,opts);
        case 'projBCD'
            loop=@(W1,H1,W2,H2) loop_projBCD(X,W1,H1,W2,H2,opts);
        case 'BCD'
            loop=@(W1,H1,W2,H2) loop_BCD(X,W1,H1,W2,H2,opts);
        otherwise
            error('Method not available')
    end

    % Main computation
    theta=opts.theta;
    if ~strcmp(opts.init,'all') 
        % Single initialization 
        [W1,H1,W2,H2]=zero_rows_pert(W1,H1,W2,H2,theta);
        [W1,H1,W2,H2,info]=loop(W1,H1,W2,H2);
        info.init=info_init;
        info.global=[];
        info.ratioinit=[];
    else
        % Multiple initialization
        [W1,H1,W2,H2]=zero_rows_pert(WH.W1_svd,WH.H1_svd,WH.W2_svd,WH.H2_svd,theta);
        [W1_svd,H1_svd,W2_svd,H2_svd,info_svd]=loop(W1,H1,W2,H2);
        
        [W1,H1,W2,H2]=zero_rows_pert(WH.W1_FS,WH.H1_FS,WH.W2_FS,WH.H2_FS,theta);
        [W1_FS,H1_FS,W2_FS,H2_FS,info_FS]=loop(W1,H1,W2,H2);
        
        [W1,H1,W2,H2]=zero_rows_pert(WH.W1_FSL,WH.H1_FSL,WH.W2_FSL,WH.H2_FSL,theta);
        [W1_FSL,H1_FSL,W2_FSL,H2_FSL,info_FSL]=loop(W1,H1,W2,H2);

        [W1,H1,W2,H2]=zero_rows_pert(WH.W1_FSR,WH.H1_FSR,WH.W2_FSR,WH.H2_FSR,theta);
        [W1_FSR,H1_FSR,W2_FSR,H2_FSR,info_FSR]=loop(W1,H1,W2,H2);

        err_vec=[info_svd.err(end),info_FS.err(end),...
            info_FSL.err(end),info_FSR.err(end)];
        method={'svd','FS','FSL','FSR'};
        W1_vec={W1_svd,W1_FS,W1_FSL,W1_FSR};
        W2_vec={W2_svd,W2_FS,W2_FSL,W2_FSR};
        H1_vec={H1_svd,H1_FS,H1_FSL,H1_FSR};
        H2_vec={H2_svd,H2_FS,H2_FSL,H2_FSR};
        info_vec={info_svd,info_FS,info_FSL,info_FSR};
        [errbest,index]=min(err_vec);
        ratioinit=zeros(4,1);
        if errbest<1e-7
            denbest=errbest;
        else
            denbest=errbest+1e-8;
        end
        for i=1:4
            ratioinit(i)=abs(errbest-err_vec(i))/denbest;
        end
        ratioinit(ratioinit<1e-4)=0;
        info=info_vec{index};
        info.init=method{index};
        info.global=info_vec;
        info.ratioinit=ratioinit;
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


