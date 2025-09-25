
    n=4;
    m=n;
    r=sqrt(n);
    eta=0;
    nmax=4;
    rng(5)
    X=randi(nmax,n,r);
    Y=randi(nmax,m,r);
    U=randi(nmax,n,r);
    V=randi(nmax,m,r);
    A=(X*Y').*(U*V')+eta*randi(nmax,n,m);
    A=A/norm(A,'fro');
    A=ones(n,m)-eye(n,m);
    %A=eye(n,m);

    block=10; %'X',X,'Y',Y,'U',U,'V',V,...
    opts=struct('rank',r,'maxit',500,'init','FS',...
        'c',1e-4,'beta',0,'Hblock',block,'Wblock',block,...
        'safestop',2,'theta',2,'tau',1e-1);
    [X_FS,Y_FS,U_FS,V_FS,err_FS]=HadDec_manifold(A,opts);
    fin_err_FS=err_FS(end);
    opts.init='Wertz';
    [X_Wertz,Y_Wertz,U_Wertz,V_Wertz,err_Wertz]=HadDec_manifold(A,opts);
    fin_err_Wertz=err_Wertz(end);
    opts.init='FS2';
    tic;
    [X_FS2,Y_FS2,U_FS2,V_FS2,err_FS2]=HadDec_manifold(A,opts);
    t_manifold=toc;
    fin_err_FS2=err_FS2(end);
    opts.init='FS3';
    [X_FS3,Y_FS3,U_FS3,V_FS3,err_FS3]=HadDec_manifold(A,opts);
    fin_err_FS3=err_FS3(end);

    close all
    lw=1.3;
    legendlabel={};
    if fin_err_FS<1e5
        semilogy(err_FS,'r-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'FS'];
    end
    if fin_err_Wertz<1e5
        semilogy(err_Wertz,'b-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'Wertz'];
    end
    if fin_err_FS2<1e5
        semilogy(err_FS2,'g-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'FS2'];
    end
    if fin_err_FS3<1e5
        semilogy(err_FS3,'m-','LineWidth',lw)
        legendlabel=[legendlabel,'FS3'];
    end
    legend(legendlabel,'Location','best')