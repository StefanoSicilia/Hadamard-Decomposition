function [W1,H1,W2,H2,info]=loop_Manopt(X,W1,H1,W2,H2,opts)
%% loop_Manopt: Main cycle for HadDec - Manopt case
% Performs Riemannian gradient descent for Hadamard decomposition. See
% HadDec_Manopt for more details.

    % Initialization and parameters
    tstart=tic;
    [m,n]=size(X);
    r=size(W1,2);
    P=struct('X1',struct('U',W1,'S',eye(r),'V',H1),'X2',...
        struct('U',W2,'S',eye(r),'V',H2));

    Mr.X1=fixedrankembeddedfactory(m,n,r);
    Mr.X2=fixedrankembeddedfactory(m,n,r);
    manifold=productmanifold(Mr);
    problem.M=manifold;
    options.verbosity=0;
    options.tolcost=opts.tol;
    options.maxiter=Inf;
    options.tolgradnorm=opts.tol;
    init_time=toc(tstart);
    options.maxtime=opts.maxtime-init_time;

    % Calling Manopt's trustregions 
    problem.cost=@(P) 0.5*norm(X-prodsvd(P.X1).*prodsvd(P.X2),'fro')^2;
    problem.egrad=@(P) grad_eval(P,X);
    problem.ehess=@(P,A) hess_eval(P,A,X,manifold);
    % warning('off', 'manopt:getHessian:approx')
    [P, ~, info] = trustregions(problem,P,options);

    % Output
    info=struct('err',sqrt(2*[info.cost]),'time',[info.time]);
    X1=P.X1;
    X2=P.X2;
    W1=X1.U*X1.S;
    H1=X1.V;
    W2=X2.U*X2.S;
    H2=X2.V;

end