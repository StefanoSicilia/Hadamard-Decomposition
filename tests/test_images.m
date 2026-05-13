%% Script to test some images 
% It takes a bit less than 4 hours. 
% The computational time is approximately given by 
% n_methods*opts.maxtime*nex*nranks

    %% Methods parameters and initializations
    opts.maxtime=240;
    opts.init='all';

    methods={'Manopt','manBCD','projBCD','BCD','TSVD'};
    n_methods=length(methods)-1;
    models={'cameraman','cat','dog1','dog2','dog3','olivettifaces'};
    nex=length(models);
    nranks=2;
    X=cell(nex,1);
    Xsvd=cell(nex,1);
    W1=cell(nex,n_methods,nranks);
    H1=cell(nex,n_methods,nranks);
    W2=cell(nex,n_methods,nranks);
    H2=cell(nex,n_methods,nranks);
    tableinfo=cell(nex,15,nranks);
    errSVD=zeros(nex,nranks);
    relerr=zeros(nex,n_methods+1,nranks);
    info=cell(nex,n_methods,nranks);
    ranks=zeros(nex,1);
    startex=1;

    %% Load the datasets
    for i=startex:nex
        switch models{i}
            case 'cameraman'
                X{i}=double(imread('cameraman.tif'));
            case 'cat'
                X{i}=double(imread('./datasets/cat.jpg'));
            case {'dog1','dog2','dog3'}
                string=['./datasets/',models{i},'.mat'];
                X{i}=load(string).data;
            case 'football'
                X{i}=load('./datasets/football.mat').A;
            case 'olivettifaces'
                X{i}=load('./datasets/olivettifaces.mat').faces;
        end
        ranks(i,2)=floor(sqrt(min(size(X{i}))));
        ranks(i,1)=round(ranks(i,2)/2);
    end

    %% Main computations
    tic
    for i=startex:nex
        fprintf('%i) %s ...',i,models{i})
        p=ranks(i,1);
        r=ranks(i,2);
        [m,n]=size(X{i});

        % HDs
        for j=1:n_methods
            fprintf(' %s ...',methods{j})
            opts.method=methods{j};
            [W1{i,j,1},H1{i,j,1},W2{i,j,1},H2{i,j,1},info{i,j,1}]=...
                HadDec(X{i},p,opts);
            relerr(i,j,1)=info{i,j,1}.err(end);
            [W1{i,j,2},H1{i,j,2},W2{i,j,2},H2{i,j,2},info{i,j,2}]=...
                HadDec(X{i},r,opts);
            relerr(i,j,2)=info{i,j,2}.err(end);
        end

        % TSVD
        p2=2*p;
        r2=min([m,n,2*r]);
        Xsvd{i}=X{i}/norm(X{i},'fro');
        [U,S,V]=svds(Xsvd{i},r2);
        errSVD(i,1)=norm(Xsvd{i}-U(:,1:p2)*S(1:p2,1:p2)*V(:,1:p2)','fro');
        errSVD(i,2)=norm(Xsvd{i}-U*S*V','fro');

        relerr(i,end,1)=errSVD(i,1);
        relerr(i,end,2)=errSVD(i,2);
        tableinfo(i,1:4,1)={models{i} m n p};
        tableinfo(i,1:4,2)={models{i} m n r};
        k=5;
        for j=1:n_methods
            tableinfo(i,k:k+1,1)={info{i,j,1}.err(end) info{i,j,1}.init};
            tableinfo(i,k:k+1,2)={info{i,j,2}.err(end) info{i,j,2}.init};
            k=k+2;
        end
        tableinfo(i,end-2,1)={errSVD(i,1)};
        tableinfo(i,end-2,2)={errSVD(i,2)};
        fprintf(' done!\n')
    end
    time_largedata=toc;


    %% Plot results for the errors
    close all
    lw=1.3;
    for j=1:n_methods
        figure
        for i=startex:nex
            rng(10*i)
            semilogy(info{i,j,2}.time,info{i,j,2}.err,'-',...
                'Color',[rand,rand,rand],'LineWidth',lw)
            hold on
        end
        title(methods{j})
        legend(models{startex:nex},'Location','best');
    end

    %% SVD improvement
    % For each dataset, find the rank for which the rank-r TSVD provides a 
    % lower relative error than the best HadDec (among manBCD and proj).

    hadbest=min(relerr(:,1:end-1,:),[],2);
    rstar_vec=2*ranks;
    ratio=zeros(nex,nranks);
    fprintf('SVD improvements checks...')
    for i=startex:nex
        for h=1:nranks
            err_star=hadbest(i,h);
            err_SVD=errSVD(i,h);
            rstar=rstar_vec(i,h);
            if err_star>err_SVD
                % rank-2r TSVD is better than rank-r HD
                while err_star>err_SVD && rstar>1
                    rstar=rstar-1;
                    [U,S,V]=svds(Xsvd{i},rstar);
                    err_SVD=norm(Xsvd{i}-U*S*V','fro');
                end
            else
                % rank-2r TSVD is worse than rank-r HD
                while err_star<err_SVD && rstar<min(size(X{i}))
                    rstar=rstar+1;
                    [U,S,V]=svds(Xsvd{i},rstar);
                    err_SVD=norm(Xsvd{i}-U*S*V','fro');
                end
            end
            rstar_vec(i,h)=rstar;
            ratio(i,h)=100*(rstar-2*ranks(i,h))/(2*ranks(i,h));
            tableinfo(i,end-1:end,h)={rstar_vec(i,h),ratio(i,h)};
        end
    end
    fprintf(' done!\n')

    %% Show images of cameraman and dog
    % Cameraman
    normX=norm(X{1},'fro');
    sq=255/normX;
    for j=1:n_methods
        figure
        imshow(sq*(W1{1,j,2}*H1{1,j,2}').*(W2{1,j,2}*H2{1,j,2}'))
        title(methods{j})
    end
    figure
    [U,S,V]=svds(Xsvd{1},2*ranks(1,2));
    imshow(sq*U*S*V'*normX)
    title(methods{end})
    figure
    imshow(sq*X{1})
    title('original')

    % Dog
    sq=255;
    A=cell(n_methods+1);
    for j=1:n_methods
        figure
        for k=1:3
            normX=norm(X{k+2},'fro');
            A{j}(:,:,k)=((W1{k+2,j,2}*H1{k+2,j,2}').*(W2{k+2,j,2}*H2{k+2,j,2}'))/normX;
        end
        imshow(sq*A{j})
        title(methods{j})
    end
    
    %%
    Xsvd_dog=zeros([size(X{3}),3]);
    figure
    for k=1:3
        [U,S,V]=svds(Xsvd{k+2},2*ranks(k+2,2));
        A{end}(:,:,k)=(U*S*V');
        normX(k)=norm(X{k+2},'fro');
        Xsvd_dog(:,:,k)=Xsvd{k+2};
    end
    imshow(sq*A{end})
    title(methods{end})
    figure
    imshow(sq*Xsvd_dog)
    title('original')

    %% Savings and percent replacement
    relerr=100*relerr;
    tableinfo(:,[5:2:end-3,end-2],:)=...
        cellfun(@(x) 100*x,tableinfo(:,[5:2:end-3,end-2],:),...
        'UniformOutput', false);

    save('./results\images_info','info')
    save('./results\images_table','tableinfo')

    for j=1:n_methods+1
        figure(j+4)
        title('')
        figstring=['./results./Images/cameraman_',methods{j},'.fig'];
        pngstring=['./results./Images/cameraman_',methods{j},'.png'];
        imwrite(frame2im(getframe(gca)),pngstring);
        title(methods{j})
        saveas(gcf,figstring)
    end
    figure(n_methods+6)
    title('')
    figstring='./results./Images/cameraman_original.fig';
    pngstring='./results./Images/cameraman_original.png';
    imwrite(frame2im(getframe(gca)),pngstring);
    title('original')
    saveas(gcf,figstring)

    for j=1:n_methods+1
        figure(j+n_methods+6)
        title('')
        figstring=['./results./Images/dog_',methods{j},'.fig'];
        pngstring=['./results./Images/dog_',methods{j},'.png'];
        imwrite(frame2im(getframe(gca)),pngstring);
        title(methods{j})
        saveas(gcf,figstring)
    end
    figure(2*n_methods+8)
    title('')
    figstring='./results./Images/dog_original.fig';
    pngstring='./results./Images/dog_original.png';
    imwrite(frame2im(getframe(gca)),pngstring);
    title('original')
    saveas(gcf,figstring)
    
