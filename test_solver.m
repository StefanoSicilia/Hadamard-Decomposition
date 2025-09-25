
    n=4;
    m=4;
    r=2;
    maxit=1;
    rng(1)
    A=eye(n,m);
    f=@(M) norm((M(1:n,1:r)*M(n+1:end,1:r)').*(M(1:n,r+1:end)*M(n+1:end,r+1:end)')-A,'fro');
    err=zeros(maxit,1);
    for j=1:maxit
        rng(j)
        X0=randi(10,n,r);
        Y0=randi(10,m,r);
        U0=randi(10,n,r);
        V0=randi(10,n,r);
        M=fminunc(f,[X0,U0;Y0,V0]);
        err(j)=f(M);
    end
    X=M(1:n,1:r);
    Y=M(n+1:end,1:r);
    U=M(1:n,r+1:end);
    V=M(n+1:end,r+1:end);
    B=(X*Y').*(U*V');
    close all
    plot(err,'r-o')