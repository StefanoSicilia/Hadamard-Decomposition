%% Script to test HadDec on synthetic data - low rank
% It takes a bit more than 3 hours

    %% Methods, parameters and structures
    m=100; 
    n=m;
    nsample=10; 
    nrank=floor(sqrt(min([m,n]))); 

    maxit=1e6;
    tol=1e-30;
    Hblock=1;
    Wblock=1;
    opts=struct('maxit',maxit,'init','all','tau',1.5,...
        'Hblock',Hblock,'Wblock',Wblock,'tol',tol);
    opts.momentum=[0.75,1,1.05,1.01,1.5];
    %opts.momentum=[0,0,0,0,1]; % algorithms without extrapolation
    opts.maxtime=40; 
    opts_1=opts; opts_1.method='manBCD';
    opts_2=opts; opts_2.method='BCD';
    opts_3=opts; opts_3.method='Manopt';
    err=zeros(nsample,nrank,4);

    % Note: total time required by HadDec is 3*opts.maxtime*nsample*nrank

    %% Apply the methods: manBCD (1), BCD (2), Manopt (3), rank-2r SVD (4)
    for i=1:nsample
        rng(i)
        for j=1:nrank
            opts_1.rank=j;
            opts_2.rank=j;
            opts_3.rank=j;
            r2=2*j;
            X=rand(m,r2)*rand(r2,n);
            % 1) manBCD
            [~,~,~,~,info_man]=HadDec(X,opts_1);
            err(i,j,1)=info_man.err(end);
            % 2) BCD
            opts.method='BCD';
            [~,~,~,~,info_BCD]=HadDec(X,opts_2);
            err(i,j,2)=info_BCD.err(end);
            % 3) Manopt
            [~,~,~,~,info_Manopt]=HadDec(X,opts_3);
            err(i,j,3)=info_Manopt.err(end);
            % 4) rank-2r SVD
            Xsvd=X/norm(X,'fro');
            [U,S,V]=svd(Xsvd);
            Sr2=S(1:r2,1:r2);
            Ur2=U(:,1:r2)*Sr2;
            Vr2=V(:,1:r2)*Sr2;
            err(i,j,4)=norm(Xsvd-Ur2*Sr2*Vr2','fro');
        end
    end 
   
    %% Store and plot results
    close all
    lw=1.3;
    err_mean=mean(err,1);
    err_std=std(err,1);
    color={'r','b','g',[0.75,0.5,0]};
    ranks=1:nrank;
    for k=1:4
        hold on;
        semilogy(ranks, err_mean(:,:,k),'-o','Color',color{k},'LineWidth',lw);
    end
    title('Mean of the errors')
    xlabel('Ranks')
    legendlabel={'manBCD','BCD','Manopt','SVD'};
    legend(legendlabel,'Location','best')



