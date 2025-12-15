%% Script to test HadDec
% It takes about 15 minutes

    %% Methods, parameters and structures
    rng(1)
    m=225;
    n=m;
    r=15;
    X=double(imread('./datasets/cameraman.jpg'));
    X=X(:,:,1);

    %% Methods parameters
    maxit=1e6;
    tol=1e-15;
    Hblock=1;
    Wblock=1;
    rng(1)
    opts=struct('rank',r,'maxit',maxit,'init','all','tau',1.5,...
        'Hblock',Hblock,'Wblock',Wblock,'tol',0.5*tol^2);
    if strcmp(opts.init,'given')
        opts.W1=W1; opts.H1=H1; opts.W2=W2; opts.H2=H2;   
    end
    opts.momentum=[0.75,1,1.05,1.01,1.5];
    %opts.momentum=[0,0,0,0,1]; % algorithms without extrapolation
    opts.maxtime=240;

    % Manifold BCD
    opts.method='manBCD';
    [W1_man,H1_man,W2_man,H2_man,info_man]=HadDec(X,opts);
    fin_err_man_BCD=info_man.err(end);

    % Wertz et al. BCD
    opts.method='BCD';
    [W1_BCD,H1_BCD,W2_BCD,H2_BCD,info_BCD]=HadDec(X,opts);
    fin_err_BCD=info_BCD.err(end);

    % Manifold Manopt
    opts.method='Manopt';
    [W1_Manopt,H1_Manopt,W2_Manopt,H2_Manopt,info_Manopt]=HadDec(X,opts);
    fin_err_Manopt=info_Manopt.err(end);

    % SVD of rank 2r
    r2=min([n,m,2*r]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svd(Xsvd);
    Sr2=S(1:r2,1:r2);
    Ur2=U(:,1:r2);
    Vr2=V(:,1:r2);
    fin_errSVDmatlab=norm(Xsvd-Ur2*Sr2*Vr2','fro');
   
    %% Display the results
    close all
    countfig=1;
    figure(countfig)
    lw=1.3;

    % Plot objective functions versus iteration 
    legendlabel={};
    if fin_err_man_BCD<1e5
        semilogy(info_man.err,'r-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'Manifold-BCD'];
    end
    if fin_err_BCD<1e5
        semilogy(info_BCD.err,'b--','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'BCD'];
    end
    if fin_err_Manopt<1e5
        semilogy(info_Manopt.err,'g-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'Manopt'];
    end   
    if fin_errSVDmatlab<1e5
        maxl=max([length(info_BCD.err),length(info_Manopt.err),...
            length(info_man.err)]);
        f=fin_errSVDmatlab;
        semilogy([0,maxl],[f,f],'-.','Color',[0.75,0.5,0],'LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'SVDmatlab'];
    end
    legend(legendlabel,'Location','best')

    % Plot objective functions versus time
    countfig=countfig+1;
    figure(countfig)
    lw=1.3;
    legendlabel={};
    if fin_err_man_BCD<1e5
        semilogy(info_man.time,info_man.err,'r-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'Manifold-BCD'];
    end
    if fin_err_Manopt<1e5
        semilogy(info_Manopt.time,info_Manopt.err,'g-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'Manopt'];
    end
    if fin_err_BCD<1e5
        semilogy(info_BCD.time,info_BCD.err,'b--','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'BCD'];
    end
    if fin_errSVDmatlab<1e5
        f=fin_errSVDmatlab;
        Mtime=max([info_BCD.time(end),info_Manopt.time(end),info_man.time(end)]);
        semilogy([0,Mtime],[f,f],'-.','Color',[0.75,0.5,0],'LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'SVDmatlab'];
    end
    legend(legendlabel,'Location','best')

    %% Show the figures if the dataset is an image
    normX=norm(X,'fro');
    sq=255/normX;
    countfig=countfig+1;
    figure(countfig)
    imshow(sq*(W1_man*H1_man').*(W2_man*H2_man'))
    title('Manifold BCD')
    countfig=countfig+1;
    figure(countfig)
    imshow(sq*(W1_BCD*H1_BCD').*(W2_BCD*H2_BCD'))
    title('BCD')
    countfig=countfig+1;
    figure(countfig)
    imshow(sq*(W1_Manopt*H1_Manopt').*(W2_Manopt*H2_Manopt'))
    title('Manopt')
    countfig=countfig+1;
    figure(countfig)
    imshow(255*Ur2*Sr2*Vr2')
    title('SVD of rank 2r')
    countfig=countfig+1;
    figure(countfig)
    imshow(sq*X)
    title('original')
    s1=ssim(X,(W1_man*H1_man').*(W2_man*H2_man'));
    sW1=ssim(X,(W1_BCD*H1_BCD').*(W2_BCD*H2_BCD'));