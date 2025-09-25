%% Script to test things on Hadamard product

    n=4;
    m=4;
    r=2;
    one=ones(r,1);
    I=eye(r);
    X=zeros(n,r);
    U=zeros(n,r);
    X(1:r^2,:)=kron(I,one);
    U(1:r^2,:)=kron(one,I);
    Y=zeros(m,r);
    V=zeros(m,r);
    Y(1:r^2,:)=kron(I,one);
    V(1:r^2,:)=kron(one,I);
    A=(X*Y').*(U*V');

    Z=column_Had(X,U);
    Q=column_Had(Y,V);
    B=Z*Q';
    err=norm(A-B,'fro');

    rng(1)
    nmax=4;
    E=randi(nmax,n,r);
    F=randi(nmax,m,r);
    C1=(E*Y').*(U*V');
    C2=(X*F').*(U*V');
    C3=(X*Y').*(E*V');
    C4=(X*Y').*(U*F');


       
    
