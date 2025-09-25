
    % Dimensions of the problem
    n=4;
    m=n;
    r=sqrt(n);
    eta=1e-5;
    nmax=4;

    % Example type
    example='randexact';
    rng(11)
    switch example
        case 'eye'
            X=eH1e(n,m);
        case 'oneseye'
            X=ones(n,m)-eye(n,m);
        case 'random'
            X=rand(n,m);
        case 'randexact'
            W1=randi(nmax,n,r);
            H1=randi(nmax,m,r);
            W2=randi(nmax,n,r);
            H2=randi(nmax,m,r);
            X=(W1*H1').*(W2*H2');
        case 'randexactpert'
            W1=randi(nmax,n,r);
            H1=randi(nmax,m,r);
            W2=randi(nmax,n,r);
            H2=randi(nmax,m,r);
            X=(W1*H1').*(W2*H2')+eta*randi(nmax,n,m);
        otherwise
            error('Example type not available.')
    end
    X=X/norm(X,'fro');

    block=10;
    opts=struct('rank',r,'maxit',50,'init','FS','tau',1,...
        'Hblock',block,'Wblock',block,'tol',1e-7);
    opts.momentum=[0,0,0,0,1];
    [W1_FS,H1_FS,W2_FS,H2_FS,err_FS]=Had_manBCD(X,opts);
    fin_err_FS=err_FS(end);
    opts.init='Wertz';
    [W1_Wertz,H1_Wertz,W2_Wertz,H2_Wertz,err_Wertz]=Had_manBCD(X,opts);
    fin_err_Wertz=err_Wertz(end);
    opts.init='FS2';
    tic;
    [W1_FS2,H1_FS2,W2_FS2,H2_FS2,err_FS2]=Had_manBCD(X,opts);
    t_manifold=toc;
    fin_err_FS2=err_FS2(end);
    opts.init='FS3';
    [W1_FS3,H1_FS3,W2_FS3,H2_FS3,err_FS3]=Had_manBCD(X,opts);
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