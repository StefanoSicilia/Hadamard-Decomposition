function [X1,X2,err]=Had_manopt4(X,opts)
%% Had_manopt4:
% Computes an approximate Hadamard decomposition X=X1.*X2, where X1 and X2
% have rank r=opts.rank. 
% It uses a power manifold Mr^2 and it optimizes through Manopt.

    [n,m]=size(X);
    r=opts.rank;
    switch opts.init
        case 'Wertz'
            M=sqrt(abs(X));
            N=sign(X).*M;
            [U1,S1,V1]=svd(M);
            W1=U1(:,1:r)*sqrt(S1(1:r,1:r));
            H1=V1(:,1:r)*sqrt(S1(1:r,1:r));
            [U2,S2,V2]=svd(N);
            W2=U2(:,1:r)*sqrt(S2(1:r,1:r));
            H2=V2(:,1:r)*sqrt(S2(1:r,1:r));
        case 'FS'
            [W1,W2,H1,H2]=Had_init(X,r);
        case 'FS2'
            [W1,W2,H1,H2]=Had_init2(X,r);
        case 'FS3'
            [W1,W2,H1,H2]=Had_init3(X,r);
        case 'given'
            W1=opts.W1;
            H1=opts.H1;
            W2=opts.W2;
            H2=opts.H2;
        otherwise
            error('Initialization not available.')
    end
    P={struct('U',W1,'S',eye(r),'V',H1),struct('U',W2,'S',eye(r),'V',H2)};


    % Manopt parameters
    tic;
    Mr=fixedrankembeddedfactory(n,m,r);
    problem.M=powermanifold(Mr,2);
    time_manifold_build=toc;

    options.verbosity=0;
    options.maxtime=opts.maxtime-time_manifold_build;

    problem.cost=@(P) 0.5*norm(X-prodsvd(P{1}).*prodsvd(P{2}),'fro')^2;
    problem.egrad=@(P) {-prodsvd(P{2}).*(X-prodsvd(P{1}).*prodsvd(P{2})),...
       -prodsvd(P{1}).*(X-prodsvd(P{1}).*prodsvd(P{2}))};
    %checkgradient(problem);
    [P, err, ~] = trustregions(problem,P,options);

    X1=P{1};
    X2=P{2};

end

function X=prodsvd(A)
    
    X=A.U*A.S*A.V';

end
