function [X1,X2,err]=Had_manopt3(X,opts)
%% Had_manopt3:
% Computes an approximate Hadamard decomposition X=X1.*X2, where X1 and X2
% have rank r=opts.rank. 
% It uses a product manifold for X1 and X2 and it optimizes through Manopt.

    [n,m]=size(X);
    r=opts.rank;
    switch opts.init
        case 'Wertz'
            manifold=sqrt(abs(X));
            N=sign(X).*manifold;
            [U1,S1,V1]=svd(manifold);
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

    % Orthogonalization of initialization
    % [Q1,R1]=qr(H1,0);
    % H1=Q1;
    % W1=W1*R1';
    % [Q2,R2]=qr(H2,0);
    % H2=Q2;
    % W2=W2*R2';


    P=struct('X1',struct('U',W1,'S',eye(r),'V',H1),'X2',...
        struct('U',W2,'S',eye(r),'V',H2));


    % Manopt parameters
    Mr.X1=fixedrankembeddedfactory(n,m,r);
    Mr.X2=fixedrankembeddedfactory(n,m,r);
    manifold=productmanifold(Mr);
    problem.M=manifold;

    options.verbosity=0;
    options.maxtime=opts.maxtime;

    problem.cost=@(P) 0.5*norm(X-prodsvd(P.X1).*prodsvd(P.X2),'fro')^2;
    problem.egrad=@(P) grad_eval(P,X);
    problem.ehess=@(P,A) hess_eval(P,A,X,manifold);
    % checkgradient(problem);
    % checkhessian(problem);
    [P, err, ~] = trustregions(problem,P,options);

    X1=P.X1;
    X2=P.X2;

end

function X=prodsvd(A)
    
    X=A.U*A.S*A.V';

end

function G=grad_eval(P,X)

    X1=prodsvd(P.X1);
    X2=prodsvd(P.X2);
    R=X-X1.*X2;
    G=struct('X1',-R.*X2,'X2',-R.*X1);

end

function H=hess_eval(P,A,X,manifold)
    
    X1=prodsvd(P.X1);
    X2=prodsvd(P.X2);
    Y=2*X1.*X2-X;
    B=manifold.tangent2ambient(P,A);
    dA1=prodsvd(B.X1);
    dA2=prodsvd(B.X2);
    H=struct('X1',dA1.*(X2.^2)+Y.*dA2,'X2',dA2.*(X1.^2)+Y.*dA1);

end


    % problem.egrad=@(P) struct('X1',-prodsvd(P.X2).*(X-prodsvd(P.X1).*prodsvd(P.X2)),...
    %     'X2',-prodsvd(P.X1).*(X-prodsvd(P.X1).*prodsvd(P.X2)));
    % problem.ehess=@(P,A) struct('X1',-(prodsvd(M.tangent2ambient(P,A).X1)).*(prodsvd(P.X2).^2)-...
    %     (X-prodsvd(P.X1).*prodsvd(P.X2)).*(prodsvd(M.tangent2ambient(P,A).X2)),...
    %     'X2',-(prodsvd(M.tangent2ambient(P,A).X2)).*(prodsvd(P.X1).^2)-...
    %     (X-prodsvd(P.X1).*prodsvd(P.X2)).*(prodsvd(M.tangent2ambient(P,A).X1)));

    % H=struct('X1',-dA1.*(X2.^2)-R.*dA2,'X2',-dA2.*(X1.^2)-R.*dA1);
    % H=struct('X1',-X2.*(X2.*dA1+X1.*dA2),'X2',-X1.*(X2.*dA1+X1.*dA2));
    % H=struct('X1',dA1.*(X2.^2)+(X1.*X2-R).*dA2,'X2',dA2.*(X1.^2)+(X1.*X2-R).*dA1);

