function [W1,H1,W2,H2,info]=loop_Manopt(X,W1,H1,W2,H2,opts)
%% loop_Manopt: Main cycle for HadDec - Manopt case
% Performs Riemannian gradient descent for Hadamard decomposition. See
% HadDec for more details.

    % Initialization and parameters
    tstart=tic;
    [m,n]=size(X);
    r1=size(W1,2);
    r2=size(W2,2);
    P=struct('X1',struct('U',W1,'S',eye(r1),'V',H1),'X2',...
        struct('U',W2,'S',eye(r2),'V',H2));

    Mr.X1=fixedrankembeddedfactory(m,n,r1);
    Mr.X2=fixedrankembeddedfactory(m,n,r2);
    manifold=productmanifold(Mr);
    problem.M=manifold;
    options.verbosity=0;
    options.tolcost=0.5*(opts.tol)^2;
    options.maxiter=Inf;
    options.tolgradnorm=0.5*(opts.tol)^2;
    init_time=toc(tstart);
    options.maxtime=opts.maxtime-init_time;

    % Calling Manopt's trustregions 
    problem.cost=@(P) 0.5*norm(X-prodsvd(P.X1).*prodsvd(P.X2),'fro')^2;
    problem.egrad=@(P) grad_eval(P,X);
    problem.ehess=@(P,A) hess_eval(P,A,X,manifold);
    % warning('off', 'manopt:getHessian:approx')
    warning('on')
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