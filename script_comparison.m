%% Script to compare Had_eBCD and Had_manBCD

    % Dimensions of the problem
    n=4;
    m=n;
    r=floor(sqrt(min(n,m)));
    eta=1e-5;
    nmax=4;
%%
    % Example type
    example='cameraman';
    rng(1)
    imageflag=0;
    switch example
        case 'eye'
            A=eye(n,m);
        case 'oneseye'
            A=ones(n,m)-eye(n,m);
        case 'random'
            A=rand(n,m);
        case 'randexact'
            X=randi(nmax,n,r);
            Y=randi(nmax,m,r);
            U=randi(nmax,n,r);
            V=randi(nmax,m,r);
            A=(X*Y').*(U*V');
        case 'randexactpert'
            X=randi(nmax,n,r);
            Y=randi(nmax,m,r);
            U=randi(nmax,n,r);
            V=randi(nmax,m,r);
            A=(X*Y').*(U*V')+eta*randi(nmax,n,m);
        case 'football'
            n=115;
            m=n;
            r=9;
            A=full(load("./datasets/football.mat").A);
        case 'cat'
            n=400;
            m=600;
            r=20;
            A=double(imread('./datasets/cat.jpg'));
            A=A+1*randn(n,m);
            imageflag=1;
        case 'cameraman'
            n=225;
            m=n;
            r=15;
            A=double(imread('./datasets/cameraman.jpg'));
            A=A(:,:,1);
            %A=A+1*randn(n,m);
            imageflag=1;
        case 'peppers'
            n=192;
            m=204;
            r=13;
            A=double(imread('./datasets/peppers.jfif'));
            A=255-A(:,:,1);
            imageflag=1;
        case 'apple'
            n=180;
            m=281;
            r=13;
            A=double(imread('./datasets/apple.jfif'));
            A=A(:,:,1);
            imageflag=1;
        case 'blender'
            n=512;
            m=n;
            r=22;
            A=double(imread('./datasets/blender.png'));
            imageflag=1;
        otherwise
            error('Example type not available.')
    end
    normA=norm(A,'fro');
    %imshow(255*A/normA);
    A=A/normA;
    %%

    maxit=30;
    tol=2e-4;
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

    % Manifold BCD
    tic;
    [X_FS,Y_FS,U_FS,V_FS,err_FS]=Had_manBCD(A,opts);
    t_FS=toc;
    fin_err_FS=err_FS(end);
    opts.init='FS2';
    tic;
    [X_FS2,Y_FS2,U_FS2,V_FS2,err_FS2]=Had_manBCD(A,opts);
    t_FS2=toc;
    fin_err_FS2=err_FS2(end);
    opts.init='FS3';
    tic;
    [X_FS3,Y_FS3,U_FS3,V_FS3,err_FS3]=Had_manBCD(A,opts);
    t_FS3=toc;
    fin_err_FS3=err_FS3(end);
    opts.init='Wertz';
    tic;
    [X_svd,Y_svd,U_svd,V_svd,err_svd]=Had_manBCD(A,opts);
    t_svd=toc;
    fin_err_svd=err_svd(end);

    % Mixed BCD
    % tic;
    % [X_FS,Y_FS,U_FS,V_FS,err_FS]=Had_BCD_mix(A,opts);
    % t_FS=toc;
    % fin_err_FS=err_FS(end);
    % opts.init='FS2';
    % tic;
    % [X_FS2,Y_FS2,U_FS2,V_FS2,err_FS2]=Had_BCD_mix(A,opts);
    % t_FS2=toc;
    % fin_err_FS2=err_FS2(end);
    % opts.init='FS3';
    % tic;
    % [X_FS3,Y_FS3,U_FS3,V_FS3,err_FS3]=Had_BCD_mix(A,opts);
    % t_FS3=toc;
    % fin_err_FS3=err_FS3(end);
    % opts.init='Wertz';
    % tic;
    % [X_svd,Y_svd,U_svd,V_svd,err_svd]=Had_BCD_mix(A,opts);
    % t_svd=toc;
    % fin_err_svd=err_svd(end);

    % Wertz BCD
    opts.init='Wertz'; opts.tol=sqrt(2*tol);
    opts.momentum=[0.75,1,1.05,1.01,1.5]; opts.maxit=maxit*(Hblock+Wblock)+2;
    %opts.momentum=[0,0,0,0,1]; opts.maxit=2*maxit*block+2;
    tic;
    [X_Wsvd,Y_Wsvd,U_Wsvd,V_Wsvd,err_Wsvd]=Had_eBCD(A,opts);
    t_Wsvd=toc;
    err_Wsvd=0.5*err_Wsvd.^2;
    fin_err_Wsvd=err_Wsvd(end);
    opts.init='FS';
    tic;
    [X_WFS,Y_WFS,U_WFS,V_WFS,err_WFS]=Had_eBCD(A,opts);
    t_WFS=toc;
    err_WFS=0.5*err_WFS.^2;
    fin_err_WFS=err_WFS(end);
    opts.init='FS2';
    tic;
    [X_WFS2,Y_WFS2,U_WFS2,V_WFS2,err_WFS2]=Had_eBCD(A,opts);
    t_WFS2=toc;
    err_WFS2=0.5*err_WFS2.^2;
    fin_err_WFS2=err_WFS2(end);
    opts.init='FS3';
    tic;
    [X_WFS3,Y_WFS3,U_WFS3,V_WFS3,err_WFS3]=Had_eBCD(A,opts);
    t_WFS3=toc;
    err_WFS3=0.5*err_WFS3.^2;
    fin_err_WFS3=err_WFS3(end);

    close all
    lw=1.3;
    legendlabel={};
    if fin_err_FS<1e5
        semilogy(err_FS,'r-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'FS'];
    end
    if fin_err_FS2<1e5
        semilogy(err_FS2,'g-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'FS2'];
    end
    if fin_err_FS3<1e5
        semilogy(err_FS3,'m-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'FS3'];
    end
    if fin_err_svd<1e5
        semilogy(err_svd,'y-','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'svd'];
    end
    if fin_err_Wsvd<1e5
        semilogy(err_Wsvd,'b--','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'Wsvd'];
    end
    if fin_err_WFS<1e5
        semilogy(err_WFS,'c--','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'WFS'];
    end
    if fin_err_WFS2<1e5
        semilogy(err_WFS2,'--','Color',[0.5,0.5,0.5],'LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'WFS2'];
    end
    if fin_err_WFS3<1e5
        semilogy(err_WFS3,'k--','LineWidth',lw)
        hold on
        legendlabel=[legendlabel,'WFS3'];
    end
    legend(legendlabel,'Location','best')

    if imageflag
        sq=255;
        figure(2)
        imshow(sq*X_svd*Y_svd')
        figure(3)
        imshow(sq*U_svd*V_svd')
        figure(4)
        imshow(sq*X_Wsvd*Y_Wsvd')
        figure(5)
        imshow(sq*U_Wsvd*V_Wsvd')
        figure(6)
        imshow(sq*(X_svd*Y_svd').*(U_svd*V_svd'))
        figure(7)
        imshow(sq*(X_Wsvd*Y_Wsvd').*(U_Wsvd*V_Wsvd'))
        figure(8)
        [U,S,V]=svd(A);
        r2=2*r;
        Sr=sqrt(S(1:r2,1:r2));
        U=U(:,1:r2)*Sr;
        V=V(:,1:r2)*Sr;
        imshow(255*U*V')
        fin_errSVDmatlab=0.5*norm(A-U*V','fro')^2;
        figure(9)
        imshow(255*A)
        s1=ssim(A,(X_svd*Y_svd').*(U_svd*V_svd'));
        sW1=ssim(A,(X_Wsvd*Y_Wsvd').*(U_Wsvd*V_Wsvd'));
    end