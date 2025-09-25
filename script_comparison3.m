%% Script to compare Had_eBCD and Had_manopt2

    % Dimensions of the problem
    n=16;
    m=n;
    r=floor(sqrt(min(n,m)));
    eta=1e-5;
    nmax=4;
%%
    % Example type
    example='randexact';
    rng(1)
    imageflag=0;
    switch example
        case 'eye'
            X=eye(n,m);
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
        case 'football'
            n=115;
            m=n;
            r=9;
            X=full(load("./datasets/football.mat").A);
        case 'cat'
            n=400;
            m=600;
            r=20;
            X=double(imread('./datasets/cat.jpg'));
            X=X+1*randn(n,m);
            imageflag=1;
        case 'cameraman'
            n=225;
            m=n;
            r=15;
            X=double(imread('./datasets/cameraman.jpg'));
            X=X(:,:,1);
            %X=X+1*randn(n,m);
            imageflag=1;
        case 'peppers'
            n=192;
            m=204;
            r=13;
            X=double(imread('./datasets/peppers.jfif'));
            X=255-X(:,:,1);
            imageflag=1;
        case 'apple'
            n=180;
            m=281;
            r=13;
            X=double(imread('./datasets/apple.jfif'));
            X=X(:,:,1);
            imageflag=1;
        case 'blender'
            n=512;
            m=n;
            r=22;
            X=double(imread('./datasets/blender.png'));
            imageflag=1;
        otherwise
            error('Example type not available.')
    end
    normX=norm(X,'fro');
    %imshow(255*X/normX);
    X=X/normX;
    %%

    maxit=20;
    tol=1e-4;
    Hblock=1;
    Wblock=1;
    rng(1)
    opts=struct('rank',r,'maxit',maxit,'init','FS','tau',2,...
        'Hblock',Hblock,'Wblock',Wblock,'tol',tol,...
        'X',randn(n,r),'Y',randn(m,r),'U',randn(n,r),'V',randn(m,r));   
        %,'X',X_Wsvd,'Y',Y_Wsvd,'U',U_Wsvd,'V',V_Wsvd);
    %opts.init='given'; 
    %opts.momentum=[0.75,1,1.05,1.01,1.5];
    opts.momentum=[0,0,0,0,1];

    % Wertz BCD
    opts.init='Wertz'; opts.tol=sqrt(2*tol);
    opts.momentum=[0.75,1,1.05,1.01,1.5]; opts.maxit=maxit*(Hblock+Wblock)+2;
    %opts.momentum=[0,0,0,0,1]; opts.maxit=2*maxit*block+2;
    tic;
    [X_Wsvd,Y_Wsvd,U_Wsvd,V_Wsvd,err_Wsvd]=Had_eBCD(X,opts);
    t_Wsvd=toc;
    err_Wsvd=0.5*err_Wsvd.^2;
    fin_err_Wsvd=err_Wsvd(end);
    opts.init='FS';
    tic;
    [X_WFS,Y_WFS,U_WFS,V_WFS,err_WFS]=Had_eBCD(X,opts);
    t_WFS=toc;
    err_WFS=0.5*err_WFS.^2;
    fin_err_WFS=err_WFS(end);
    opts.init='FS2';
    tic;
    [X_WFS2,Y_WFS2,U_WFS2,V_WFS2,err_WFS2]=Had_eBCD(X,opts);
    t_WFS2=toc;
    err_WFS2=0.5*err_WFS2.^2;
    fin_err_WFS2=err_WFS2(end);
    opts.init='FS3';
    tic;
    [X_WFS3,Y_WFS3,U_WFS3,V_WFS3,err_WFS3]=Had_eBCD(X,opts);
    t_WFS3=toc;
    err_WFS3=0.5*err_WFS3.^2;
    fin_err_WFS3=err_WFS3(end);

    % Manopt approach Had_manopt2
    warning('off', 'manopt:getHessian:approx')
    opts.init='Wertz'; opts.tol=tol; opts.maxit=maxit*(Hblock+Wblock)+2;
    %opts.momentum=[0,0,0,0,1]; opts.maxit=2*maxit*block+2;
    tic;
    [X1_Msvd,X2_Msvd,err_Msvd]=Had_manopt2(X,opts);
    t_Msvd=toc;
    fin_err_Msvd=err_Msvd(end);
    opts.init='FS';
    tic;
    [X1_MFS,X2_MFS,err_MFS]=Had_manopt2(X,opts);
    t_MFS=toc;
    fin_err_MFS=err_MFS(end);
    opts.init='FS2';
    tic;
    [X1_MFS2,X2_MFS2,err_MFS2]=Had_manopt2(X,opts);
    t_MFS2=toc;
    fin_err_MFS2=err_MFS2(end);
    opts.init='FS3';
    tic;
    [X1_MFS3,X2_MFS3,err_MFS3]=Had_manopt2(X,opts);
    t_MFS3=toc;
    fin_err_MFS3=err_MFS3(end);

    close all
    lw=1.3;
    legendlabel={};
    if fin_err_Wsvd<1e5
        semilogy(err_Wsvd,'y-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'Wsvd'];
    end
    if fin_err_WFS<1e5
        semilogy(err_WFS,'r-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'WFS'];
    end
    if fin_err_WFS2<1e5
        semilogy(err_WFS2,'g-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'WFS2'];
    end
    if fin_err_WFS3<1e5
        semilogy(err_WFS3,'m-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'WFS3'];
    end
    if fin_err_Msvd<1e5
        semilogy(err_Msvd,'b--','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'Msvd'];
    end
    if fin_err_MFS<1e5
        semilogy(err_MFS,'c--','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'MFS'];
    end
    if fin_err_MFS2<1e5
        semilogy(err_MFS2,'--','Color',[0.5,0.5,0.5],'LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'MFS2'];
    end
    if fin_err_MFS3<1e5
        semilogy(err_MFS3,'k--','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'MFS3'];
    end
    legend(legendlabel,'Location','best')

    if imageflag
        sq=255;
        figure(2)
        imshow(sq*(X_svd*Y_svd').*(U_svd*V_svd'))
        figure(3)
        imshow(sq*prodsvd(X1_Msvd).*prodsvd(X2_Msvd))
        figure(4)
        [W2,S,H2]=svd(X);
        r2=2*r;
        Sr=sqrt(S(1:r2,1:r2));
        W2=W2(:,1:r2)*Sr;
        H2=H2(:,1:r2)*Sr;
        imshow(255*W2*H2')
        fin_errSVDmatlab=0.5*norm(X-W2*H2','fro')^2;
        figure(5)
        imshow(255*X)
        s1=ssim(X,(X_svd*Y_svd').*(U_svd*V_svd'));
        sW1=ssim(X,prodsvd(X1_Msvd).*prodsvd(X2_Msvd));
    end


    function X=prodsvd(A)
    
        X=A.U*A.S*A.V';

    end