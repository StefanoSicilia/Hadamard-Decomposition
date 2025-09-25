function [W,H]=Had_Dec(A,r,opts)

    n=size(A,1);
    m=size(A,2);
    maxit=opts.maxit;
    tol=opts.tol;
    h=opts.h;
    rng(1)
    W=rand(n,r)*rand(r,m);
    H=rand(n,r)*rand(r,m);

    j=1;
    M=A-W.*H;
    normA=norm(A,'fro');
    err=norm(M,'fro')/normA;
    while j<maxit && err>=tol
        W=W-h*Projrank(M.*H,W);
        H=H-h*Projrank(M.*W,H);
        M=A-W.*H;
        err=norm(M,'fro')/normA
        j=j+1;
    end
end