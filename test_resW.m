%% Testing resW and resgradW

    n=4;
    m=5;
    r=2;
    eta=0;
    nmax=4;
    rng(1)
    X=randi(nmax,n,r);
    Y=randi(nmax,m,r);
    U=randi(nmax,n,r);
    V=randi(nmax,m,r);
    A=(X*Y').*(U*V')+eta*randi(nmax,n,m);

    RW={cell(n,1),cell(n,1)};
    RH={cell(m,1),cell(m,1)};
    for i=1:n
        RW{1}{i}=X(i,:);
        RW{2}{i}=U(i,:);
    end
    for i=1:m
        RH{1}{i}=Y(i,:);
        RH{2}{i}=V(i,:);
    end
    W=face_split(X,U);
    H=face_split(Y,V);
    a=resW(RW,A,H);
    b=resW(RH,A',W);
    normGa=norm(resgradW(RW,A,H),'fro');
    normGb=norm(resgradW(RH,A',W),'fro');

    
    
    