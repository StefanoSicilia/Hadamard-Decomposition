%% Script to test HadDec

    %% Choice of the example type
    m=4;
    n=4;
    r=floor(sqrt(min(m,n)));
    eta=1e-4;
    nmax=4;
    example='random';
    rng(1)
    switch example
        case 'random'
            X=rand(m,n);
        case 'randexact'
            W1=randi(nmax,m,r);
            H1=randi(nmax,n,r);
            W2=randi(nmax,m,r);
            H2=randi(nmax,n,r);
            X=(W1*H1').*(W2*H2');
        case 'randexactpert'
            W1=randi(nmax,m,r);
            H1=randi(nmax,n,r);
            W2=randi(nmax,m,r);
            H2=randi(nmax,n,r);
            X=(W1*H1').*(W2*H2')+eta*randi(nmax,m,n);
        case 'zeroperts'
            X=zeros(m,n);
            X(m,:)=ones(n,1);
            W1=[ones(1,r); zeros(m-1,r)];
            H1=[ones(1,r); zeros(n-1,r)];
            W2=[ones(1,r); zeros(m-1,r)];
            H2=[ones(1,r); zeros(n-1,r)];
        otherwise
            error('Example type not available.')
    end

    %% Methods parameters
    maxit=1e6;
    tol=1e-16;
    Iter_W=4;
    Iter_H=4;
    opts=struct('maxit',maxit,'init','all','tau',0.95,'theta',1e-4,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',0,'noloops',1);
    if strcmp(opts.init,'given')
        opts.W1=W1; opts.H1=H1; opts.W2=W2; opts.H2=H2;   
    end
    opts.momentum=[0.75,1,1.05,1.01,1.5,0.65];
    opts.maxtime=0.5;
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
    close all
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
