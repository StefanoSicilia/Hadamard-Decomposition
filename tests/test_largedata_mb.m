%% Script to test some large datasets - Manopt and BCD
% It takes a bit less than 6 hours.
% The computational time is approximately given by 
% (runcases{1}+runcases{2})*opts.maxtime

    %% Methods parameters and initializations
    opts.maxtime=800;
    opts.init='all';

    methods={'manBCD','projBCD'};
    n_methods=length(methods);
    models={'NG20','classic','la12','ohscal','sports','la1','la2',...
        'reviews','k1b','hitech','tr41','tr45','tr11','tr23'};
    nex=length(models);
    W1=cell(nex,n_methods);
    H1=cell(nex,n_methods);
    W2=cell(nex,n_methods);
    H2=cell(nex,n_methods);
    err=cell(nex,10);
    errSVD=zeros(nex,1);
    relerr=zeros(nex,n_methods);
    info=cell(nex,n_methods);
    runcases={2:14,5:14};
    startex=min([runcases{1},runcases{2}]);

    %% Dimensions of the datasets
    % dims_m=[19949 7094 11162 2340 2301 4069 8580 3204 6279 3075 414 ...
    %     204 878 690];
    % dims_n=[43586 41681 11465 21839 10080 18483 14870 31472 31472 31472 ...
    %     6429 5832 7454 8261];

    %% Main computations
    tstart=tic;
    ranks=[20 4 6 10 7 6 6 5 6 6 10 10 9 6];
   
    % Manopt
    for i=startex:nex
        fprintf('%i) %s ...',i,models{i})
        r=ranks(i);
        string=['./datasets/',models{i},'.mat'];
        X=load(string).dtm;
        [m,n]=size(X);
        for j=1:n_methods
            if ismember(i,runcases{j})
                fprintf(' %s ...',methods{j})
                opts.method=methods{j};
                [W1{i,j},H1{i,j},W2{i,j},H2{i,j},info{i,j}]=HadDec(X,r,opts);
                relerr(i,j)=info{i,j}.err(end);
            else
                relerr(i,j)=NaN;
                info{i,j}.err=NaN;
                info{i,j}.init='none';
                info{i,j}.time=0;
            end
        end
        r2=min([n,m,2*r]);
        Xsvd=X/norm(X,'fro');
        [U,S,V]=svds(Xsvd,r2);
        errSVD(i)=norm(Xsvd-U*S*V','fro');
        err(i,1:4)={m n r errSVD(i)};
        k=5;
        for j=1:n_methods
            err(i,k:k+1)={info{i,j}.err(end) info{i,j}.init};
            k=k+2;
        end
        fprintf(' done!\n')
    end
    time_largedata=toc(tstart);


    %% Plot results for the errors
    close all
    lw=1.3;
    for j=1:n_methods
        figure
        for i=startex:nex
            rng(10*i)
            semilogy(info{i,j}.time,info{i,j}.err,'-',...
                'Color',[rand,rand,rand],'LineWidth',lw)
            hold on
        end
        title(methods{j})
        legend(models{runcases{j}},'Location','best');
    end

    %% SVD improvement
    % For each dataset, find the rank for which the rank-r TSVD provides a 
    % lower relative error than the best HadDec (among manBCD and proj).

    hadbest=min(relerr,[],2);
    rstar_vec=2*ranks;
    ratio=zeros(nex,1);
    for i=startex:nex
        err_star=hadbest(i);
        err_SVD=errSVD(i);
        rstar=rstar_vec(i);
        string=['./datasets/',models{i},'.mat'];
        X=load(string).dtm;
        Xsvd=X/norm(X,'fro');
        if err_star>err_SVD
            % rank-2r TSVD is better than rank-r HD
            while err_star>err_SVD
                rstar=rstar-1;
                [U,S,V]=svds(Xsvd,rstar);
                err_SVD=norm(Xsvd-U*S*V','fro');
            end
        else
            % rank-2r TSVD is worse than rank-r HD
            while err_star<err_SVD
                rstar=rstar+1;
                [U,S,V]=svds(Xsvd,rstar);
                err_SVD=norm(Xsvd-U*S*V','fro');
            end
        end
        rstar_vec(i)=rstar;
        ratio(i)=100*(rstar-2*ranks(i))/(2*ranks(i));
        err(i,9:10)={rstar_vec(i),ratio(i)};
    end

    %% Savings
    % save('./results\largedata_mb_info','info')
    % save('./results\largedata_mb_err','err')