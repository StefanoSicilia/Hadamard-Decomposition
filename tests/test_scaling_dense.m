%% Script to test performances of HD with different sizes
% It takes about 8 hours

    %% Methods, parameters and structures
    n_dim=7;
    n_sample=10; 
    incr=6;
    
    init_vec={'SVD-based','FS','best'};
    maxdim_vec=[13 13 13 12];
    methods={'Manopt','manBCD','projBCD','BCD','TSVD'};
    n_init=length(init_vec);
    n_methods=length(methods)-1;
    opts.maxtime=200;
    opts.toler=1e-4;
    opts.tol=1e-3;
    start_dim=1;

    W1=cell(n_sample,n_methods,n_dim,n_init);
    W2=cell(n_sample,n_methods,n_dim,n_init);
    H1=cell(n_sample,n_methods,n_dim,n_init);
    H2=cell(n_sample,n_methods,n_dim,n_init);
    relerr=zeros(n_sample,n_methods+1,n_dim,n_init);
    times=zeros(n_sample,n_methods+1,n_dim,n_init);
    info=cell(n_sample,n_methods,n_dim,n_init);

    X=cell(n_sample,n_dim);
    normX=zeros(n_sample,n_dim);
    Xsvd=cell(n_sample,n_dim);
    U=cell(n_sample,n_dim);
    S=cell(n_sample,n_dim);
    V=cell(n_sample,n_dim);
    Ur2=cell(n_sample,n_dim);
    Sr2=cell(n_sample,n_dim);
    Vr2=cell(n_sample,n_dim);

    %% Apply the methods 
    tstart=tic;
    for i=1:n_sample
        fprintf('Element %i of the sample \n',i)
        for k=start_dim:n_dim
            rng(i*(n_dim+1)+k)
            exponent=k+incr;
            m=2^exponent;
            n=m;
            r=exponent;
            r2=2*r;
            X{i,k}=(rand(m,r)*rand(r,n)).*(rand(m,r)*rand(r,n));

            % TSVD 
            normX(i,k)=norm(X{i,k},'fro');
            Xsvd{i,k}=X{i,k}/normX(i,k);
            tic;
            [U{i,k},S{i,k},V{i,k}]=svd(Xsvd{i,k});
            times(i,end,k,:)=toc;
            Sr2{i,k}=S{i,k}(1:r2,1:r2);
            Ur2{i,k}=U{i,k}(:,1:r2);
            Vr2{i,k}=V{i,k}(:,1:r2);
            relerr(i,end,k,:)=...
                norm(Xsvd{i,k}-Ur2{i,k}*Sr2{i,k}*Vr2{i,k}','fro');

            for h=1:n_init-1

                opts.init=init_vec{h};

                % Manopt
                j=1;
                opts.method=methods{j};
                opts.toler=0;
                [W1{i,j,k,h},H1{i,j,k,h},W2{i,j,k,h},H2{i,j,k,h},info{i,j,k,h}]=...
                HadDec(X{i,k},r,opts);
                relerr(i,j,k,h)=info{i,j,k,h}.err(end);
                times(i,j,k,h)=info{i,j,k,h}.time(end);
    
                % manBCD and projBCD
                for j=2:3
                    opts.method=methods{j};
                    opts.toler=1e-4;
                    [W1{i,j,k,h},H1{i,j,k,h},W2{i,j,k,h},H2{i,j,k,h},info{i,j,k,h}]=...
                    HadDec(X{i,k},r,opts);
                    relerr(i,j,k,h)=info{i,j,k,h}.err(end);
                    times(i,j,k,h)=info{i,j,k,h}.time(end);
                end
                
                % BCD
                j=4;
                if exponent<=maxdim_vec(j)
                    opts.method=methods{j};
                    opts.toler=1e-4;
                    [W1{i,j,k,h},H1{i,j,k,h},W2{i,j,k,h},H2{i,j,k,h},info{i,j,k,h}]=...
                    HadDec(X{i,k},r,opts);
                    relerr(i,j,k,h)=info{i,j,k,h}.err(end);
                    times(i,j,k,h)=info{i,j,k,h}.time(end);
                else
                    relerr(i,j,k,h)=NaN;
                    times(i,j,k,h)=NaN;
                    info{i,j,k,h}=struct('err',[],'time',[],'init',[],'global',[],'ratioinit',[]);
                    W1{i,j,k,h}=[];
                    H1{i,j,k,h}=[];
                    W2{i,j,k,h}=[];
                    H2{i,j,k,h}=[];
                end
            end

            h=3;
            for j=1:n_methods
                [relerr(i,j,k,h),index]=min(relerr(i,j,k,1:2));
                times(i,j,k,h)=times(i,j,k,index);
                info{i,j,k,h}=info{i,j,k,index};
                W1{i,j,k,h}=W1{i,j,k,index};
                H1{i,j,k,h}=H1{i,j,k,index};
                W2{i,j,k,h}=W2{i,j,k,index};
                H2{i,j,k,h}=H2{i,j,k,index};
            end
            
        end
    end 
    t_global=toc(tstart); 
 
    %% Store and plot results
    close all
    lw=1.3;
    err_mean=squeeze(mean(relerr,1));
    err_std=squeeze(std(relerr,1));
    times_mean=squeeze(mean(times,1));
    times_std=squeeze(std(times,1));
    plotsettings={'magenta','red','blue','black',[0.75,0.5,0]};
    linestyles={'-','--','-.',':','--'};
    linemarks={'o','d','*','s','^'};
    xlabels=(1:n_dim)+incr;
    h=3;
    figure
    for j=1:n_methods+1
        hold on;
        semilogy(xlabels, err_mean(j,h),linemarks{j},'Color',plotsettings{j},...
            'LineWidth',lw,'LineStyle',linestyles{j});
    end
    title('Mean of the errors')
    xlabel('$k$','Interpreter','latex')
    xticks(xlabels)
    legend(methods,'Location','best')

    figure
    for j=1:n_methods+1
        hold on;
        semilogy(xlabels, times_mean(j,h),linemarks{j},'Color',plotsettings{j},...
            'LineWidth',lw,'LineStyle',linestyles{j});
    end
    title('Mean of the times')
    xlabel('$k$','Interpreter','latex')
    xticks(xlabels)
    legend(methods,'Location','best')

    %% Save
    % save('./final_results\scaling_dense_info','info')
    % save('./final_results\scaling_dense_err','relerr')
    % save('./final_results\scaling_dense_times','times')
    
    