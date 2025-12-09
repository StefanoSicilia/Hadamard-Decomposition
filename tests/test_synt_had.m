%% Script to test HadDec on synthetic data - low rank Hadamard
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
    res=zeros(nsample,nrank,5);

    % Note: total time required by HadDec is 3*opts.maxtime*nsample*nrank

    %% Apply the methods: manBCD (1), BCD (2), Manopt (3), rank-2r SVD (4)
    for i=1:nsample
        rng(i)
        for j=1:nrank
            opts_1.rank=j;
            opts_2.rank=j;
            opts_3.rank=j;
            r2=2*j;
            X=(rand(m,j)*rand(j,n)).*(rand(m,j)*rand(j,n));
            % 1) manBCD
            [~,~,~,~,info_man]=HadDec(X,opts_1);
            res(i,j,1)=info_man.err(end);
            % 2) BCD
            opts.method='BCD';
            [~,~,~,~,info_BCD]=HadDec(X,opts_2);
            res(i,j,2)=info_BCD.err(end);
            % 3) Manopt
            [~,~,~,~,info_Manopt]=HadDec(X,opts_3);
            res(i,j,3)=info_Manopt.err(end);
            % 4) rank-2r SVD
            Xsvd=X/norm(X,'fro');
            [U,S,V]=svd(Xsvd,'econ');
            Sr2=S(1:r2,1:r2);
            Ur2=U(:,1:r2);
            Vr2=V(:,1:r2);
            res(i,j,4)=norm(Xsvd-Ur2*Sr2*Vr2','fro');
            % 5) Svd matches Had error
            errhad=min(res(i,j,1:3));
            errsvd=res(i,j,4);
            k=r2;
            while errsvd>=errhad && k<m
                k=k+1;
                [U,S,V]=svd(Xsvd,'econ');
                Sk=S(1:k,1:k);
                Uk=U(:,1:k);
                Vk=V(:,1:k);
                errsvd=norm(Xsvd-Uk*Sk*Vk','fro');
            end
            res(i,j,5)=k;
        end
    end 
   
    %% Store and plot results
    close all
    lw=1.3;
    res_mean=mean(res,1);
    res_std=std(res,1);
    color={'r','b','g',[0.75,0.5,0]};
    ranks=1:nrank;
    for k=1:4
        hold on;
        semilogy(ranks, res_mean(:,:,k),'-o','Color',color{k},'LineWidth',lw);
    end
    title('Mean of the errors')
    xlabel('Ranks')
    legendlabel={'manBCD','BCD','Manopt','SVD'};
    legend(legendlabel,'Location','best')



