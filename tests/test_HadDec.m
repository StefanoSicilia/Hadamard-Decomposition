%% Script to test HadDec

    %% Choice of the example type
    m=16;
    n=m;
    r=floor(sqrt(min(m,n)));
    r1=r;
    r2=2*r-r1;
    
    nmax=4;
    example='random';
    rng(1)
    switch example
        case 'random'
            X=rand(m,n);
        case 'randexact'
            W1=randi(nmax,m,r1);
            H1=randi(nmax,n,r1);
            W2=randi(nmax,m,r2);
            H2=randi(nmax,n,r2);
            X=(W1*H1').*(W2*H2');
        case 'randexactpert'
            eta=1e-4;
            W1=randi(nmax,m,r1);
            H1=randi(nmax,n,r1);
            W2=randi(nmax,m,r2);
            H2=randi(nmax,n,r2);
            X=(W1*H1').*(W2*H2')+eta*randi(nmax,m,n);
        case 'zeroperts'
            X=zeros(m,n);
            X(m,:)=ones(n,1);
            W1=[ones(1,r1); zeros(m-1,r1)];
            H1=[ones(1,r1); zeros(n-1,r1)];
            W2=[ones(1,r2); zeros(m-1,r2)];
            H2=[ones(1,r2); zeros(n-1,r2)];
        otherwise
            error('Example type not available.')
    end

    %% Methods parameters
    maxit=1e6;
    tol=1e-16;
    Iter_W=2;
    Iter_H=2;
    opts=struct('maxit',maxit,'init','all','tau',0.95,'theta',1e-4,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',0,...
        'noloops',1,'rescale',1);
    if strcmp(opts.init,'given')
        opts.W1=W1; opts.H1=H1; opts.W2=W2; opts.H2=H2;   
    end
    opts.rescale=0;
    opts.momentum=[0.75,1,1.05,1.01,1.5,0.6];
    opts.maxtime=1;
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
        [W1{k},H1{k},W2{k},H2{k},info{k}]=HadDec(X,[r1 r2],opts);
        times(k)=toc;
        relerr(k)=info{k}.err(end);
    end

    %% TSVD of rank 2r
    rsvd=min([n,m,r1+r2]);
    Xsvd=X/norm(X,'fro');
    tic;
    [U,S,V]=svd(Xsvd);
    times(end)=toc;
    Srsvd=S(1:rsvd,1:rsvd);
    Ursvd=U(:,1:rsvd);
    Vrsvd=V(:,1:rsvd);
    relerr(end)=norm(Xsvd-Ursvd*Srsvd*Vrsvd','fro');
   
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
