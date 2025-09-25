function [X1,X2,err]=Had_manopt2(X,opts)
%% Had_manopt2:
% Computes an approximate Hadamard decomposition X=X1.*X2, where X1 and X2
% have rank r=opts.rank. 
% It uses a 2 block coordinate descend algorithm for X1 and X2 and 
% optimizes each rank-r matrix through Manopt.

    [n,m]=size(X);
    r=opts.rank;
    maxit=opts.maxit;
    tol=opts.tol;
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
    % X1=struct('L',W1,'R',H1);
    % X2=struct('L',W2,'R',H2);
    X1=struct('U',W1,'S',eye(r),'V',H1);
    X2=struct('U',W2,'S',eye(r),'V',H2);


    % Manopt parameters
    %Mr=fixedrankfactory_2factors(n,m,r);
    Mr=fixedrankembeddedfactory(n,m,r);
    problem1.M=Mr;
    problem2.M=Mr;
    err=zeros(maxit,1);
    err(1)=0.5*norm(X-prodsvd(X1).*prodsvd(X2),'fro')^2;
    normX=norm(X,'fro');
    a=0.5*norm(X-(W1*H1').*(W2*H2'),'fro')^2/normX;
    b=0.5*norm(X-(prodsvd(X1)).*(prodsvd(X2)),'fro')^2/normX;

    options.verbosity=0;
    j=2;
    currerr=1+tol;
    
    while j<=maxit && currerr>=tol

        % X1-update
        % problem1.cost=@(X1) 0.5*norm(X-(X1.L*X1.R').*(X2.L*X2.R'),'fro')^2;
        % problem1.egrad=@(X1) -(X2.L*X2.R').*(X-(X1.L*X1.R').*(X2.L*X2.R'));
        % problem1.ehess=@(X1,A) (X2.L*X2.R').*(X2.L*X2.R')*A;
        problem1.cost=@(X1) 0.5*norm(X-prodsvd(X1).*prodsvd(X2),'fro')^2;
        problem1.egrad=@(X1) -prodsvd(X2).*(X-prodsvd(X1).*prodsvd(X2));
        %problem1.ehess=@(X1,A) (prodsvd(X2).^2).*problem1.M.tangent2ambient(X1,A);
        %checkhessian(problem1);
        [X1, rescost, ~] = trustregions(problem1,X1,options);
        err(j)=rescost;
        j=j+1;

        % X2-update
        % problem2.cost=@(X2) 0.5*norm(X-(X1.L*X1.R').*(X2.L*X2.R'),'fro')^2;
        % problem2.egrad=@(X2) -(X1.L*X1.R').*(X-(X1.L*X1.R').*(X2.L*X2.R'));
        % problem2.ehess=@(X2,A) (X1.L*X1.R').*(X1.L*X1.R')*A;
        problem2.cost=@(X2) 0.5*norm(X-prodsvd(X1).*prodsvd(X2),'fro')^2;
        problem2.egrad=@(X2) -prodsvd(X1).*(X-prodsvd(X1).*prodsvd(X2));
        %problem1.ehess=@(X2,A) A.*(prodsvd(X1).^2);
        %checkgradient(problem1);
        [X2,rescost, ~] = trustregions(problem2,X2,options);
        err(j)=rescost;
        j=j+1;

        currerr=err(j-1);

    end
    err=err(1:j-1);

end

function X=prodsvd(A)
    
    X=A.U*A.S*A.V';

end
