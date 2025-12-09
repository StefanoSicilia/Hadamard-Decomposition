%% Script to test some large datasets 
% It takes ~3 hours

    %% Methods parameters and initializations
    maxit=1e6;
    tol=1e-30;
    Hblock=1;
    Wblock=1;
    opts=struct('maxit',maxit,'init','sparse','tau',1.5,...
        'Hblock',Hblock,'Wblock',Wblock,'tol',tol,'method','manBCDsparse');
    opts.momentum=[0.75,1,1.05,1.01,1.5];
    %opts.momentum=[0,0,0,0,1]; % algorithms without extrapolation
    opts.maxtime=800;

    nex=14;
    err=cell(nex,6);
    info=cell(nex,1);
    i=1;
    legendlabel={};

    % Note: total time required by HadDec is nex*opts.maxtime

    %% NG20
    fprintf('%i) NG20 ...',i)
    m=19949;
    n=43586;
    opts.rank=20;
    X=load("./datasets/NG20.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'NG20'];
    i=i+1;

    %% classic
    fprintf('%i) classic ...',i)
    m=7094;
    n=41681;
    opts.rank=4;
    X=load("./datasets/classic.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'classic'];
    i=i+1;

    %% ohscal
    fprintf('%i) ohscal ...',i)
    m=11162;
    n=11465;
    opts.rank=10;
    X=load("./datasets/ohscal.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'ohscal'];
    i=i+1;

    %% k1b
    fprintf('%i) k1b ...',i)
    m=2340;
    n=21839;
    opts.rank=6;
    X=load("./datasets/k1b.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'k1b'];
    i=i+1;

    %% Hitech
    fprintf('%i) Hitech ...',i)
    m=2301;
    n=10080;
    opts.rank=6;
    X=load("./datasets/hitech.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'hitech'];
    i=i+1;

    %% Reviews
    fprintf('%i) Reviews ...',i)
    m=4069;
    n=18483;
    opts.rank=5;
    X=load("./datasets/reviews.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'reviews'];
    i=i+1;

    %% Sports
    fprintf('%i) Sports ...',i)
    m=8580;
    n=14870;
    opts.rank=7;
    X=load("./datasets/sports.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'sports'];
    i=i+1;

    %% la1
    fprintf('%i) la1 ...',i) 
    m=3204;
    n=31472;
    opts.rank=6;
    X=load("./datasets/la1.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'la1'];
    i=i+1;

    %% la12
    fprintf('%i) la12 ...',i)
    m=6279;
    n=31472;
    opts.rank=6;
    X=load("./datasets/la12.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'la12'];
    i=i+1;

    %% la2
    fprintf('%i) la2 ...',i)
    m=3075;
    n=31472;
    opts.rank=6;
    X=load("./datasets/la2.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'la2'];
    i=i+1;

    %% tr11
    fprintf('%i) tr11 ...',i)
    m=414;
    n=6429;
    opts.rank=9;
    X=load("./datasets/tr11.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'tr11'];
    i=i+1;

    %% tr23
    fprintf('%i) tr23 ...',i)
    m=204;
    n=5832;
    opts.rank=6;
    X=load("./datasets/tr23.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'tr23'];
    i=i+1;

    %% tr41
    fprintf('%i) tr41 ...',i)
    m=878;
    n=7454;
    opts.rank=10;
    X=load("./datasets/tr41.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'tr41'];
    i=i+1;

    %% tr45
    fprintf('%i) tr45 ...',i)
    m=690;
    n=8261;
    opts.rank=10;
    X=load("./datasets/tr45.mat").dtm;
    [~,~,~,~,info{i}]=HadDec(X,opts);
    r2=min([n,m,2*opts.rank]);
    Xsvd=X/norm(X,'fro');
    [U,S,V]=svds(Xsvd,r2);
    errSVD=norm(Xsvd-U*S*V','fro');
    err(i,:)={m n opts.rank errSVD info{i}.err(end) info{i}.init};
    fprintf(' done!\n')
    legendlabel=[legendlabel,'tr45'];

    %% Plot results for the errors
    close all
    figure(1)
    lw=1.3;
    for i=1:nex
        rng(10*i)
        semilogy(info{i}.time,info{i}.err,'-','Color',[rand,rand,rand],...
            'LineWidth',lw)
        hold on
    end
    legend(legendlabel,'Location','best');

    %% Ratio
    ratio=zeros(nex,1);
    for i=1:nex
        ratio(i)=(info{i}.err(end-1)-info{i}.err(end))/info{i}.err(end);
    end

    %% SVD improvement
    % Use the following lines for each dataset to found the rank for which
    % SVD provides a lower relative error than HadDec.

    % X=load("./datasets/classic.mat").dtm;
    % r2=23;
    % Xsvd=X/norm(X,'fro');
    % [U,S,V]=svds(Xsvd,r2);
    % errSVD=norm(Xsvd-U*S*V','fro')

    % Ranks found to improve HadDec relative error:
    % [22 13 33 25 26 18 26 25 22 22 24 15 36 28]




