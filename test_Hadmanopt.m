%% Not working :(
    n=4;
    m=4;
    r=2;
    eta=1e-5;
    nmax=4;
    rng(1)
    X=randi(nmax,n,r);
    Y=randi(nmax,m,r);
    U=randi(nmax,n,r);
    V=randi(nmax,m,r);
    A=(X*Y').*(U*V')+eta*randi(nmax,n,m);

    opts=struct('rank',r,'maxit',10,'init','FS');
    [X,U,Y,V,err]=Had_manopt(A,opts);