
    n=4;
    m=4;
    r=2;
    eta=0;
    nmax=4;
    rng(1)
    X=randi(nmax,n,r);
    Y=randi(nmax,m,r);
    U=randi(nmax,n,r);
    V=randi(nmax,m,r);
    A=(X*Y').*(U*V')+eta*randi(nmax,n,m);
    A=ones(n,m)-eye(n,m);
    %A=eye(n,m);

    opts=struct('rank',r,'maxit',1000,'init','FS3','h',0.1,...
        'alpha',5e-3);
    [X_FS,Y_FS,U_FS,V_FS,err_FS]=Had_riemopt(A,opts);
    fin_err_FS=err_FS(end);
    opts.init='Wertz';
    [X_Wertz,Y_Wertz,U_Wertz,V_Wertz,err_Wertz]=Had_riemopt(A,opts);
    fin_err_Wertz=err_Wertz(end);
    opts.init='FS2';
    [X_FS2,Y_FS2,U_FS2,V_FS2,err_FS2]=Had_riemopt(A,opts);
    fin_err_FS2=err_FS2(end);
    opts.init='FS3';
    [X_FS3,Y_FS3,U_FS3,V_FS3,err_FS3]=Had_riemopt(A,opts);
    fin_err_FS3=err_FS3(end);

    close all
    lw=1.3;
    semilogy(err_FS,'r-','LineWidth',lw)
    hold on
    semilogy(err_Wertz,'b-','LineWidth',lw)
    hold on
    semilogy(err_FS2,'g-','LineWidth',lw)
    hold on
    semilogy(err_FS3,'m-','LineWidth',lw)
    legend('FS','Wertz','FS2','FS3','Location','best')