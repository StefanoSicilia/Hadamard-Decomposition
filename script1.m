%% Script to test Had_Dec

    n=10;
    m=10;
    r=8;
    rng(1)
    A=rand(n,m);
    opts=struct('maxit',10000,'h',1,'tol',1e-6);
    [W,H]=Had_Dec(A,r,opts);
    norm(A-W.*H,'fro')/norm(A,'fro')