%% Script to extrapolation and rescaling on HD for synthetic data
% Tests HD on three different types of random synthetic data (general,
% low-rank and Hadamard-decomposable) with an without extrapolation and
% with or without rescaling (see options 'momentum' and 'rescale' in
% HadDec).
% It takes approximately 2 hours.

    %% Methods, parameters and structures
    m=400; 
    n=m;
    nsample=10;
    
    opts.init='all';
    methods={'Manopt','manBCD','projBCD','BCD','TSVD'};
    types={'general-rank','low-rank','Had-low-rank'};
    kinds={'none','extrapolation','rescale','both'};
    extra_vec_mp=[0,0,0,0,1; 0.25,1,1.05,1.01,1.5;...
        0,0,0,0,1;0.25,1,1.05,1.01,1.5];
    extra_vec_bcd=[0,0,0,0,1; 0.75,1,1.05,1.01,1.5;...
        0,0,0,0,1;0.75,1,1.05,1.01,1.5];
    scale_vec=[0,0,1,1];
    ranks=[10 15 20];
    nmethods=length(methods)-1;
    ntype=length(types);
    nrank=length(ranks);
    nkind=length(kinds);
    maxtimes=[10 10 10]; 

    W1=cell(nsample,nrank,nmethods,ntype,nkind);
    W2=cell(nsample,nrank,nmethods,ntype,nkind);
    H1=cell(nsample,nrank,nmethods,ntype,nkind);
    H2=cell(nsample,nrank,nmethods,ntype,nkind);
    relerr=zeros(nsample,nrank,nmethods+1,ntype,nkind);
    info=cell(nsample,nrank,nmethods+1,ntype,nkind);

    X=zeros(m,n,nsample,nrank,ntype);
    normX=zeros(nsample,nrank,ntype);
    Xsvd=zeros(m,n,nsample,nrank,ntype);
    U=cell(nsample,nrank,ntype);
    S=cell(nsample,nrank,ntype);
    V=cell(nsample,nrank,ntype);
    Ursvd=cell(nsample,nrank,ntype);
    Srsvd=cell(nsample,nrank,ntype);
    Vrsvd=cell(nsample,nrank,ntype);

    %% Apply the methods 
    for i=1:nsample
        fprintf('Element %i of the sample \n',i)
        for k=1:nrank
            rng(i*(nrank+1)+k)
            r=ranks(k);
            rsvd=2*r;
            X(:,:,i,k,1)=rand(m,n);
            X(:,:,i,k,2)=rand(m,rsvd)*rand(rsvd,n);
            X(:,:,i,k,3)=(rand(m,r)*rand(r,n)).*(rand(m,r)*rand(r,n));

            for h=1:ntype

                opts.maxtime=maxtimes(h);
                j=1; 
                l=1;
                [W1{i,k,j,h,l},H1{i,k,j,h,l},W2{i,k,j,h,l},...
                    H2{i,k,j,h,l},info{i,k,j,h,l}]=...
                    HadDec(X(:,:,i,k,h),r,opts);
                relerr(i,k,j,h,l)=info{i,k,j,h,l}.err(end);

                % HDs
                for j=2:3
                    opts.method=methods{j};
                    for l=1:nkind
                        opts.momentum=extra_vec_mp(l,:);
                        opts.rescale=scale_vec(l);
                        [W1{i,k,j,h,l},H1{i,k,j,h,l},W2{i,k,j,h,l},...
                            H2{i,k,j,h,l},info{i,k,j,h,l}]=...
                            HadDec(X(:,:,i,k,h),r,opts);
                        relerr(i,k,j,h,l)=info{i,k,j,h,l}.err(end);
                    end
                end

                j=4;
                opts.method=methods{j};
                for l=[1,2]
                    opts.momentum=extra_vec_bcd(l,:);
                    [W1{i,k,j,h,l},H1{i,k,j,h,l},W2{i,k,j,h,l},...
                        H2{i,k,j,h,l},info{i,k,j,h,l}]=...
                        HadDec(X(:,:,i,k,h),r,opts);
                    relerr(i,k,j,h,l)=info{i,k,j,h,l}.err(end);
                end

                 
                % TSVD 
                normX(i,k,h)=norm(X(:,:,i,k,h),'fro');
                Xsvd(:,:,i,k,h)=X(:,:,i,k,h)/normX(i,k,h);
                [U{i,k,h},S{i,k,h},V{i,k,h}]=svd(Xsvd(:,:,i,k,h));
                Srsvd{i,k,h}=S{i,k,h}(1:rsvd,1:rsvd);
                Ursvd{i,k,h}=U{i,k,h}(:,1:rsvd);
                Vrsvd{i,k,h}=V{i,k,h}(:,1:rsvd);
                relerr(i,k,end,h,1)=norm(Xsvd(:,:,i,k,h)-...
                    Ursvd{i,k,h}*Srsvd{i,k,h}*Vrsvd{i,k,h}','fro');

            end
        end
    end 
   
    %% Store and plot results
    err_mean=squeeze(mean(relerr,1))*100;
    err_std=squeeze(std(relerr,1))*100;

    %% Save
    % save('./results\extrscale_info','info')
    % save('./results\extrscale_err','relerr')
    
    