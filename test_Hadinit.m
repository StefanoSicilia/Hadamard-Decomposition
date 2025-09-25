%% Testing Had_init

    n=16;
    m=20;
    r=4;
    maxit=1;
    err_new=zeros(maxit,1);
    err_WVG=zeros(maxit,1);
    epsilon=1e-2;

    for j=1:maxit
        rng(j)
        Xtrue=randi(10,n,r);
        Ytrue=randi(10,m,r);
        Utrue=randi(10,n,r);
        Vtrue=randi(10,m,r);
        E=randi(10,n,m);
        A=(Xtrue*Ytrue').*(Utrue*Vtrue')+epsilon*E/norm(E,'fro');
        [X,U,Y,V]=Had_init(A,r);
        err_new(j)=norm(A-face_split(X,U)*face_split(Y,V)','fro')/norm(A,'fro');
   
        M=sqrt(abs(A));
        N=sign(A).*M;
        [U1,S1,V1]=svd(M);
        W1=U1(:,1:r)*sqrt(S1(1:r,1:r));
        H1=V1(:,1:r)*sqrt(S1(1:r,1:r));
        [U2,S2,V2]=svd(N);
        W2=U2(:,1:r)*sqrt(S2(1:r,1:r));
        H2=V2(:,1:r)*sqrt(S2(1:r,1:r));
        err_WVG(j)=norm(A-(W1*H1').*(W2*H2'),'fro')/norm(A,'fro');
    end

    close all
    plot(err_new,'r-o')
    hold on
    plot(err_WVG,'b-o');
    legend('New approach','Wertz et al.')
    score=sum(err_WVG>err_new)/maxit;