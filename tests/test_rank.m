%% Script to test HadDec on low-rank synthetic data for rank decrease
% It takes about 2.5 hours

    %% Methods, parameters and structures
    m=1000; 
    n=m;
    nsample=10; 
    nrank=30;
      
    maxit=1e6;
    tol=1e-16;
    Iter_W=2;
    Iter_H=2;
    opts=struct('maxit',maxit,'init','all','tau',0.95,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',0,'noloop',1);
    opts.momentum=[0.75,1,1.05,1.01,1.5];
    opts.maxtime=30;
    opts.method='Manopt';
    err=zeros(nsample,nrank);
    init=cell(nsample,nrank);    

    %% Apply the Manopt method
    for i=1:nsample
        rng(i)
        for j=1:nrank
            r2=2*j;
            X=rand(m,r2)*rand(r2,n);
            [~,~,~,~,info_Manopt]=HadDec(X,j,opts);
            err(i,j)=info_Manopt.err(end);
            init{i,j}=info_Manopt.init;
        end
    end 
   
    %% Store and plot results
    lw=1.3;
    err_mean=mean(err,1);
    err_std=std(err,1);
    color={'g',[0.75,0.5,0]};
    ranks=1:nrank;
    figure
    semilogy(ranks, err_mean,'-o','Color',color{1},'LineWidth',lw);
    title('Mean of the errors')
    xlabel('Ranks')
    legendlabel={'Manopt'};
    legend(legendlabel,'Location','best')



