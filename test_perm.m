%% Test for some permutations

    r=2;
    n=r^2;
    m=n;
    nmax=4;
    rng(1)
    rng(21)
    X=[randi(nmax,r,r);eye(r)];
    U=[ones(r);randi(nmax,r,r)];
    Y=[randi(nmax,r,r);eye(r)];
    V=[eye(r);randi(nmax,r,r)];
    A=(X*Y').*(U*V');

    P=[0,0,0,1;0,0,1,0;0,1,0,0;1,0,0,0];
    Q=[0,0,1,0;1,0,0,0;0,1,0,0;0,0,0,1];
    % norm(Q*face_split(X,U)-face_split(Q*X,Q*U),'fro')
    norm(face_split(Y,V)'*Q'-face_split(Q*Y,Q*V)','fro')
    norm(P*A*Q'-face_split(P*X,P*U)*face_split(Q*Y,Q*V)')