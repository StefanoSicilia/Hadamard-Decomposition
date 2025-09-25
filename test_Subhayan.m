

    n=4;
    m=n;
    r=1;
    nmax=4;
    for i=1:1071
        rng(i)
        M=randi(nmax,n,r);
        A=M*randi(nmax,r,m);
        B=randi(nmax,n,r)*M';
        a=max(vecnorm(A,2,1));
        b=max(vecnorm(B,2,1));
        ab=max(vecnorm(A*B,2,1));
        if ab==sqrt(n)*a*b
            disp('aaa')
            i
        end
    end
    