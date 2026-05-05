%% Script to test error trend of HD for synthetic data
% It takes a bit more than 8 hours.
% The computational time is approximately given by 
% sum(maxtimes)*n_methods*nsample*ntype

    %% Methods, parameters and structures
    m=400; 
    n=m;
    n_sample=10; 
    
    maxit=1e6;
    tol=1e-16;
    Iter_W=2;
    Iter_H=2;
    opts=struct('maxit',maxit,'init','all','tau',0.95,'theta',1e-4,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',1,'noloops',1);
    opts.momentum=[0.75,1,1.05,1.01,1.5,0.6];
    methods={'Manopt','manBCD','projBCD','BCD','TSVD'};
    types={'general-rank','low-rank','Had-low-rank'};
    ranks=[10 15 20];
    n_methods=length(methods)-1;
    ntype=length(types);
    nrank=length(ranks);
    maxtimes=[40 100 100];

    W1=cell(n_sample,nrank,n_methods,ntype);
    W2=cell(n_sample,nrank,n_methods,ntype);
    H1=cell(n_sample,nrank,n_methods,ntype);
    H2=cell(n_sample,nrank,n_methods,ntype);
    relerr=zeros(n_sample,nrank,n_methods+1,ntype);
    info=cell(n_sample,nrank,n_methods+1,ntype);
    best_init=zeros(n_sample,nrank,n_methods,4,ntype);

    X=zeros(m,n,n_sample,nrank,ntype);
    normX=zeros(n_sample,nrank,ntype);
    Xsvd=zeros(m,n,n_sample,nrank,ntype);
    U=cell(n_sample,nrank,ntype);
    S=cell(n_sample,nrank,ntype);
    V=cell(n_sample,nrank,ntype);
    Ur2=cell(n_sample,nrank,ntype);
    Sr2=cell(n_sample,nrank,ntype);
    Vr2=cell(n_sample,nrank,ntype);

    %% Apply the methods 
    tstart=tic;
    for i=1:n_sample
        fprintf('Element %i of the sample \n',i)
        for k=1:nrank
            rng(i*(nrank+1)+k)
            r=ranks(k);
            r2=2*r;
            X(:,:,i,k,1)=rand(m,n);
            X(:,:,i,k,2)=rand(m,r2)*rand(r2,n);
            X(:,:,i,k,3)=(rand(m,r)*rand(r,n)).*(rand(m,r)*rand(r,n));

            for h=1:ntype

                opts.maxtime=maxtimes(h);

                % HDs
                for j=1:n_methods
                    opts.method=methods{j};
                    [W1{i,k,j,h},H1{i,k,j,h},W2{i,k,j,h},H2{i,k,j,h},info{i,k,j,h}]=...
                    HadDec(X(:,:,i,k,h),r,opts);
                    relerr(i,k,j,h)=info{i,k,j,h}.err(end);
                    best_init(i,k,j,:,h)=double(info{i,k,j,h}.ratioinit==0);
                end
                 
                % TSVD 
                normX(i,k,h)=norm(X(:,:,i,k,h),'fro');
                Xsvd(:,:,i,k,h)=X(:,:,i,k,h)/normX(i,k,h);
                [U{i,k,h},S{i,k,h},V{i,k,h}]=svd(Xsvd(:,:,i,k,h));
                Sr2{i,k,h}=S{i,k,h}(1:r2,1:r2);
                Ur2{i,k,h}=U{i,k,h}(:,1:r2);
                Vr2{i,k,h}=V{i,k,h}(:,1:r2);
                relerr(i,k,end,h)=...
                    norm(Xsvd(:,:,i,k,h)-Ur2{i,k,h}*Sr2{i,k,h}*Vr2{i,k,h}','fro');

                % Comparison between HD and TSVD
                hadbest=min(relerr(i,k,1:end-1,h));
                rstar=r2;
                err_star=hadbest;
                err_SVD=relerr(i,k,end,h);
                if err_star>err_SVD
                    % rank-2r TSVD is better than rank-r HD
                    while err_star>err_SVD && rstar>1
                        rstar=rstar-1;
                        err_SVD=norm(Xsvd(:,:,i,k,h)-U{i,k,h}(:,1:rstar)*S{i,k,h}(1:rstar,1:rstar)*V{i,k,h}(:,1:rstar)','fro');
                    end
                else
                    % rank-2r TSVD is worse than rank-r HD
                    while err_star<err_SVD && rstar<rank(Xsvd(:,:,i,k,h))
                        rstar=rstar+1;
                        err_SVD=norm(Xsvd(:,:,i,k,h)-U{i,k,h}(:,1:rstar)*S{i,k,h}(1:rstar,1:rstar)*V{i,k,h}(:,1:rstar)','fro');
                    end
                end
                ratio=100*(rstar-r2)/(r2);
                info{i,k,end,h}=[err_SVD,rstar,ratio];
            end
        end
    end 
    t_global=toc(tstart);

    % Best initializations
    init_vec=permute(squeeze(sum(best_init)),[3 1 2 4]);
   
    %% Store and plot results
    close all
    lw=1.3;
    err_mean=squeeze(mean(relerr,1));
    err_std=squeeze(std(relerr,1));
    plotsettings={'magenta','red','blue','black',[0.75,0.5,0]};
    for h=1:ntype
        figure
        for j=1:n_methods+1
            hold on;
            semilogy(ranks, err_mean(:,j,h),'-o','Color',plotsettings{j},'LineWidth',lw);
        end
        title(['Mean of the errors ',types{h},' case'])
        xlabel('Ranks')
        legend(methods,'Location','best')
    end

    %% Save
    % save('./results\trend_info','info')
    % save('./results\trend_err','relerr')
    % save('./results\trend_init','init_vec')
    
    