%% Script to test extrapolation and rescaling for HD
% The computational time is approximately given by 
% sum(maxtimes)*n_methods*nsample*ntype

    %% Methods, parameters and structures
    m=100; 
    n=m;
    r=10;
    rsvd=2*r;
    n_sample=5; 
    
    maxit=1e6;
    tol=1e-16;
    Iter_W=2;
    Iter_H=2;
    opts=struct('maxit',maxit,'init','all','tau',0.95,'theta',1e-4,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',0,'noloops',1);
    momentums={[0,0,0,0,1,0.6],[0.75,1,1.05,1.01,1.5,0.6],...
        [0,0,0,0,1,0.6],[0.75,1,1.05,1.01,1.5,0.6]};
    rescales=[0 0 1 1];
    methods={'Manopt','manBCD','projBCD','BCD','TSVD'};
    types={'none','extrapolation','rescaling','both'};
    n_methods=length(methods)-1;
    ntype=length(types);
    maxtimes=5*[1 1 1 1];

    W1=cell(n_sample,n_methods,ntype);
    W2=cell(n_sample,n_methods,ntype);
    H1=cell(n_sample,n_methods,ntype);
    H2=cell(n_sample,n_methods,ntype);
    relerr=zeros(n_sample,n_methods+1,ntype);
    info=cell(n_sample,n_methods+1,ntype);

    X=zeros(m,n,n_sample);
    normX=zeros(n_sample,1);
    Xsvd=zeros(m,n,n_sample);
    U=cell(n_sample,1);
    S=cell(n_sample,1);
    V=cell(n_sample,1);
    Ur2=cell(n_sample,1);
    Sr2=cell(n_sample,1);
    Vr2=cell(n_sample,1);

    %% Apply the methods 
    tstart=tic;
    for i=1:n_sample
        fprintf('Element %i of the sample \n',i)
        rng(i)
        X(:,:,i)=rand(m,n);

        for h=1:ntype

            opts.maxtime=maxtimes(h);
            opts.momentum=momentums{h};
            opts.rescale=rescales(h);

            % HDs
            for j=1:n_methods
                fprintf('%s ...',methods{j})
                opts.method=methods{j};
                [W1{i,j,h},H1{i,j,h},W2{i,j,h},H2{i,j,h},info{i,j,h}]=...
                HadDec(X(:,:,i),r,opts);
                relerr(i,j,h)=info{i,j,h}.err(end);
            end
            fprintf('\n')
        end

        % TSVD 
        normX(i)=norm(X(:,:,i),'fro');
        Xsvd(:,:,i)=X(:,:,i)/normX(i);
        [U{i},S{i},V{i}]=svd(Xsvd(:,:,i));
        Sr2{i}=S{i}(1:rsvd,1:rsvd);
        Ur2{i}=U{i}(:,1:rsvd);
        Vr2{i}=V{i}(:,1:rsvd);
        relerr(i,end,:)=norm(Xsvd(:,:,i)-Ur2{i}*Sr2{i}*Vr2{i}','fro')*ones(ntype,1);

    end 
    t_global=toc(tstart);
   
    %% Mean and std of the results
    err_mean=squeeze(mean(relerr,1));
    err_std=squeeze(std(relerr,1));

    
    