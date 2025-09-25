%% Testing manopt things

    n=4;
    m=4;
    r=2;
    Mr.X1=fixedrankembeddedfactory(n,m,r);
    Mr.X2=fixedrankembeddedfactory(n,m,r);
    M=productmanifold(Mr);
    problem.M=M;
    X=M.rand();
    Y=M.rand();
    U=M.proj(X,struct('X1',rand(n,r)*rand(r,m),'X2',rand(n,r)*rand(r,m)));
    V=M.tangent2ambient(X,U)