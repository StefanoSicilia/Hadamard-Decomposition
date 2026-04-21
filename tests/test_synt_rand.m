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
    opts=struct('maxit',maxit,'init','all','tau',0.95,'theta',1e-4,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',0,'noloops',1);
    opts.momentum=[0.75,1,1.05,1.01,1.5,0.6];
    opts.maxtime=10; 
    methods={'Manopt','manBCD','projBCD','BCD','TSVD'};
    n_methods=length(methods)-1;
    W1=cell(n_sample,nrank,n_methods);
    W2=cell(n_sample,nrank,n_methods);
    H1=cell(n_sample,nrank,n_methods);
    H2=cell(n_sample,nrank,n_methods);
    times=zeros(n_sample,nrank,n_methods);
    relerr=zeros(n_sample,nrank,n_methods+1);
    info=cell(n_sample,nrank,n_methods+1);
    best_init=zeros(n_sample,nrank,n_methods,4);

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

        % HD
        for j=1:n_methods
            opts.method=methods{j};
            for r=1:nrank
                tic;
                [W1{i,r,j},H1{i,r,j},W2{i,r,j},H2{i,r,j},info{i,r,j}]=...
                HadDec(X(:,:,i),r+1,opts);
                times(i,r,j)=toc;
                relerr(i,r,j)=info{i,r,j}.err(end);
                best_init(i,r,j,:)=double(info{i,r,j}.ratioinit==0);
            end
        end

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

            % Comparison between HD and TSVD
            hadbest=min(relerr(i,r,1:end-1));
            rstar=r2;
            err_star=hadbest;
            err_SVD=relerr(i,r,end);
            if err_star>err_SVD
                % rank-2r TSVD is better than rank-r HD
                while err_star>err_SVD && rstar>1
                    rstar=rstar-1;
                    err_SVD=norm(Xsvd(:,:,i)-U{i}(:,1:rstar)*S{i}(1:rstar,1:rstar)*V{i}(:,1:rstar)','fro');
                end
            else
                % rank-2r TSVD is worse than rank-r HD
                while err_star<err_SVD && rstar<min(m,n)
                    rstar=rstar+1;
                    err_SVD=norm(Xsvd(:,:,i)-U{i}(:,1:rstar)*S{i}(1:rstar,1:rstar)*V{i}(:,1:rstar)','fro');
                end
            end
            ratio=100*(rstar-r2)/(r2);
            info{i,r,end}=[err_SVD,rstar,ratio];
        end
    end 

    % Best initializations
    init_vec=permute(squeeze(sum(best_init)),[3 1 2]);
   
    %% Store and plot results
    close all
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

    %% Save the results (optional)
    % save('./results\synt_rand_info','info')
    % save('./results\synt_rand_err','relerr')
    % save('./results\synt_rand_init','init_vec')




