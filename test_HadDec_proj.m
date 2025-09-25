
    n=4;
    m=4;
    r=2;
    eta=0;
    nmax=4;
    rng(4)
    X=randi(nmax,n,r);
    Y=randi(nmax,m,r);
    U=randi(nmax,n,r);
    V=randi(nmax,m,r);
    A=(X*Y').*(U*V')+eta*randi(nmax,n,m);
    % A=ones(n,m)-eye(n,m);
    % A=eye(n,m);

    W=face_split(X,U);
    H=face_split(Y,V);
    opts=struct('rank',r,'maxit',100,'init','FS3');%,'given','W',W,'H',H);
    [W_FS,H_FS,err_FS]=HadDec_proj(A,opts);
    fin_err_FS=err_FS(end);
    opts.init='Wertz';
    [W_Wertz,H_Wertz,err_Wertz]=HadDec_proj(A,opts);
    fin_err_Wertz=err_Wertz(end);
    opts.init='FS2';
    tic;
    [W_FS2,H_FS2,err_FS2]=HadDec_proj(A,opts);
    t_manifold=toc;
    fin_err_FS2=err_FS2(end);
    opts.init='FS3';
    [W_FS3,H_FS3,err_FS3]=HadDec_proj(A,opts);
    fin_err_FS3=err_FS3(end);

    close all
    lw=1.3;
    semilogy(err_FS,'r-','LineWidth',lw)
    hold on
    %semilogy(err_Wertz,'b-','LineWidth',lw)
    hold on
    %semilogy(err_FS2,'g-','LineWidth',lw)
    hold on
    %semilogy(err_FS3,'m-','LineWidth',lw)
    legend('FS','Wertz','FS2','FS3','Location','best')