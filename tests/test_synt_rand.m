%% Script to test HadDec on synthetic data - random data
% It takes a bit more than 1 hour.
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
    opts=struct('maxit',maxit,'init','all','tau',0.95,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',0,'noloop',1);
    opts.momentum=[0.75,1,1.05,1.01,1.5];
    opts.maxtime=10; 
    methods={'Manopt','manBCD','proj','BCD','TSVD'};
    n_methods=length(methods)-1;
    W1=cell(n_sample,nrank,n_methods);
    W2=cell(n_sample,nrank,n_methods);
    H1=cell(n_sample,nrank,n_methods);
    H2=cell(n_sample,nrank,n_methods);
    times=zeros(n_sample,nrank,n_methods);
    relerr=zeros(n_sample,nrank,n_methods+1);
    info=cell(n_sample,nrank,n_methods+1);

    X=zeros(m,n,n_sample);
    normX=zeros(n_sample,1);
    Xsvd=zeros(m,n,n_sample);
    U=cell(n_sample,nrank);
    S=cell(n_sample,nrank);
    V=cell(n_sample,nrank);
    Ur2=cell(n_sample,nrank);
    Sr2=cell(n_sample,nrank);
    Vr2=cell(n_sample,nrank);

    %% Apply the methods 
    for i=1:n_sample
        rng(i)
        X(:,:,i)=rand(m,n);

        % TSVD
        normX(i)=norm(X(:,:,i),'fro');
        Xsvd(:,:,i)=X(:,:,i)/normX(i);
        [U{i},S{i},V{i}]=svd(Xsvd(:,:,i));
        for r=1:nrank
            r2=2*(r+1);
            Sr2{i,r}=S{i}(1:r2,1:r2);
            Ur2{i,r}=U{i}(:,1:r2);
            Vr2{i,r}=V{i}(:,1:r2);
            relerr(i,r,end)=norm(Xsvd(:,:,i)-Ur2{i,r}*Sr2{i,r}*Vr2{i,r}','fro');
        end

        % HD
        for j=1:n_methods
            opts.method=methods{j};
            for r=1:nrank
                tic;
                [W1{i,r,j},H1{i,r,j},W2{i,r,j},H2{i,r,j},info{i,r,j}]=...
                HadDec(X(:,:,i),r+1,opts);
                times(i,r,j)=toc;
                relerr(i,r,j)=info{i,r,j}.err(end);
            end
        end
    end 

   
    %% Store and plot results
    figure
    lw=1.3;
    err_mean=mean(relerr,1);
    err_std=std(relerr,1);
    plotsettings={'magenta','red','blue','black',[0.75,0.5,0]};
    ranks=2:nrank+1;
    for j=1:n_methods+1
        hold on;
        semilogy(ranks, err_mean(:,:,j),'-o','Color',plotsettings{j},'LineWidth',lw);
    end
    title('Mean of the errors')
    xlabel('Ranks')
    legendlabel=methods;
    legend(legendlabel,'Location','best')




