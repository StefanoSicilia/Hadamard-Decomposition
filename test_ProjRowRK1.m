%% Testing ProjRowRK1

    n=4;
    m=4;
    r=2;
    A=ones(n,m); A(1,1)=2;
    M=ProjRowRK1(A);