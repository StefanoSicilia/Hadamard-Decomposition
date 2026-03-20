%% HadDec on cameraman image
% It takes almost 20 minutes.
% The computational time is approximately given by 
% n_methods*opts.maxtime

    %% Methods, parameters and structures
    m=225;
    n=m;
    r=15; 
    X=double(imread('./datasets/cameraman.jpg'));
    X=X(:,:,1);

   %% Methods parameters
    maxit=1e6;
    tol=1e-16;
    Iter_W=3;
    Iter_H=3;
    rng(1)
    opts=struct('maxit',maxit,'init','all','tau',0.95,'theta',1e-4,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',1,'noloops',1);
    if strcmp(opts.init,'given')
        opts.W1=W1; opts.H1=H1; opts.W2=W2; opts.H2=H2;   
    end
    opts.momentum=[0.75,1,1.05,1.01,1.5,0.6];
    opts.maxtime=240;
    methods={'Manopt','manBCD','projBCD','BCD','TSVD'};
    n_methods=length(methods)-1;
    W1=cell(n_methods,1);
    W2=cell(n_methods,1);
    H1=cell(n_methods,1);
    H2=cell(n_methods,1);
    times=zeros(n_methods+1,1);
    relerr=zeros(n_methods+1,1);
    info=cell(n_methods+1,1);

    %% Rank-r HDs
    for k=1:n_methods
        opts.method=methods{k};
        tic;
        [W1{k},H1{k},W2{k},H2{k},info{k}]=HadDec(X,r,opts);
        times(k)=toc;
        relerr(k)=info{k}.err(end);
    end

    %% TSVD of rank 2r
    r2=min([n,m,2*r]);
    Xsvd=X/norm(X,'fro');
    tic;
    [U,S,V]=svd(Xsvd);
    times(end)=toc;
    Sr2=S(1:r2,1:r2);
    Ur2=U(:,1:r2);
    Vr2=V(:,1:r2);
    relerr(end)=norm(Xsvd-Ur2*Sr2*Vr2','fro');
   
    %% Display the results
    %close all
    figure
    lw=1.3;
    plotsettings={'m-','r-','b-','k-'};

    % Plot objective functions versus iteration 
    for k=1:n_methods
        semilogy(info{k}.err,plotsettings{k},'LineWidth',lw)
        hold on 
    end
    f=relerr(end);
    maxl=1;
    for k=1:n_methods
        err_length=length(info{k}.err);
        if err_length>maxl
            maxl=err_length;
        end
    end
    semilogy([0,maxl],[f,f],'-.','Color',[0.75,0.5,0],'LineWidth',lw)
    hold on
    legend(methods,'Location','best')

    % Plot objective functions versus time
    figure
    for k=1:n_methods
        semilogy(info{k}.time,info{k}.err,plotsettings{k},'LineWidth',lw)
        hold on
    end
    f=relerr(end);
    maxt=0;
    for k=1:n_methods
        time_method=info{k}.time(end);
        if time_method>maxt
            maxt=time_method;
        end
    end
    semilogy([0,maxt],[f,f],'-.','Color',[0.75,0.5,0],'LineWidth',lw)
    hold on
    legend(methods,'Location','best')

    %% Show the figures of the image approximations
    normX=norm(X,'fro');
    sq=255/normX;
    for k=1:n_methods
        figure
        imshow(sq*(W1{k}*H1{k}').*(W2{k}*H2{k}'))
        title(methods{k})
    end
    figure
    imshow(sq*Ur2*Sr2*Vr2'*normX)
    title(methods{end})
    figure
    imshow(sq*X)
    title('original')

