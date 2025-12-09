%% Script to test HadDec on low-rank synthetic data for rank decrease
% It takes about 2.5 hours

    %% Methods, parameters and structures
    m=1000; 
    n=m;
    nsample=10; 
    nrank=30;
      
    maxit=1e6;
    tol=1e-30;
    Hblock=1;
    Wblock=1;
    opts=struct('maxit',maxit,'init','all','tau',1.5,...
        'Hblock',Hblock,'Wblock',Wblock,'tol',tol);
    opts.momentum=[0.75,1,1.05,1.01,1.5];
    %opts.momentum=[0,0,0,0,1]; % algorithms without extrapolation
    opts.maxtime=30;
    opts.method='Manopt';
    err=zeros(nsample,nrank);
    init=cell(nsample,nrank);    

    % Note: total time required by HadDec is opts.maxtime*nsample*nrank 

    %% Apply the Manopt method
    for i=1:nsample
        rng(i)
        for j=1:nrank
            opts.rank=j;
            r2=2*j;
            X=rand(m,r2)*rand(r2,n);
            [~,~,~,~,info_Manopt]=HadDec(X,opts);
            err(i,j)=info_Manopt.err(end);
            init{i,j}=info_Manopt.init;
        end
    end 
   
    %% Store and plot results
    close all
    lw=1.3;
    err_mean=mean(err,1);
    err_std=std(err,1);
    color={'g',[0.75,0.5,0]};
    ranks=1:nrank;
    semilogy(ranks, err_mean(:,:,k),'-o','Color',color{k},'LineWidth',lw);
    title('Mean of the errors')
    xlabel('Ranks')
    legendlabel={'Manopt'};
    legend(legendlabel,'Location','best')



