%% Script to test HadDec on synthetic data - low-rank in Hadamard sense
% It takes a bit more than 4 hours.
% The computational time is approximately given by 
% n_methods*opts.maxtime*nsample*nrank

    %% Methods, parameters and structures
    m=100; 
    n=m;
    n_sample=10; 
    nrank=floor(sqrt(min([m,n])))-1;    
    
    maxit=1e6;
    tol=1e-16;
    Iter_W=2;
    Iter_H=2;
    opts=struct('maxit',maxit,'init','all','tau',0.95,'theta',1e-4,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',0,'noloops',1);
    opts.momentum=[0.75,1,1.05,1.01,1.5,0.6];
    opts.maxtime=40; 
    methods={'Manopt','manBCD','projBCD','BCD','TSVD'};
    n_methods=length(methods)-1;
    W1=cell(n_sample,nrank,n_methods);
    W2=cell(n_sample,nrank,n_methods);
    H1=cell(n_sample,nrank,n_methods);
    H2=cell(n_sample,nrank,n_methods);
    times=zeros(n_sample,nrank,n_methods);
    relerr=zeros(n_sample,nrank,n_methods+1);
    info=cell(n_sample,nrank,n_methods+1);

    X=zeros(m,n,n_sample,nrank);
    normX=zeros(n_sample,nrank);
    Xsvd=zeros(m,n,n_sample,nrank);
    U=cell(n_sample,nrank);
    S=cell(n_sample,nrank);
    V=cell(n_sample,nrank);
    Ur2=cell(n_sample,nrank);
    Sr2=cell(n_sample,nrank);
    Vr2=cell(n_sample,nrank);

    %% Apply the methods 
    for i=1:n_sample
        for r=1:nrank
            rng(i*(nrank+1)+r)
            r2=2*(r+1);
            X(:,:,i,r)=(rand(m,r+1)*rand(r+1,n)).*(rand(m,r+1)*rand(r+1,n));

            % TSVD
            normX(i,r)=norm(X(:,:,i,r),'fro');
            Xsvd(:,:,i,r)=X(:,:,i,r)/normX(i,r);
            [U{i,r},S{i,r},V{i,r}]=svd(Xsvd(:,:,i,r));
            Sr2{i,r}=S{i,r}(1:r2,1:r2);
            Ur2{i,r}=U{i,r}(:,1:r2);
            Vr2{i,r}=V{i,r}(:,1:r2);
            relerr(i,r,end)=norm(Xsvd(:,:,i,r)-Ur2{i,r}*Sr2{i,r}*Vr2{i,r}','fro');
    
            % HD
            for j=1:n_methods
                opts.method=methods{j};
                tic;
                [W1{i,r,j},H1{i,r,j},W2{i,r,j},H2{i,r,j},info{i,r,j}]=...
                HadDec(X(:,:,i,r),r+1,opts);
                times(i,r,j)=toc;
                relerr(i,r,j)=info{i,r,j}.err(end);
            end
        end
    end 

   
    %% Store and plot results
    figure
    lw=1.3;
    err_mean=squeeze(mean(relerr,1));
    err_std=squeeze(std(relerr,1));
    plotsettings={'magenta','red','blue','black',[0.75,0.5,0]};
    ranks=2:nrank+1;
    for j=1:n_methods+1
        hold on;
        semilogy(ranks, err_mean(:,j),'-o','Color',plotsettings{j},'LineWidth',lw);
    end
    title('Mean of the errors')
    xlabel('Ranks')
    legendlabel=methods;
    legend(legendlabel,'Location','best')



