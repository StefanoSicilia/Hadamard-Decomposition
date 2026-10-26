%% Script to test HadDec on simple examples for rank-r HD
% To run the test, follow the steps and choose among the possible options.

    %% 1) Choose the example
    % Possible choices: 'rand-gen','rand-lr','rand-hlr','rand-sparse',
    % 'rand-sparse-hlr' or 'cameraman'
    example='rand-hlr';
    
    rng(1)
    switch example
        case 'rand-gen' % Random generic example
            m=9;
            n=m;
            r=floor(sqrt(min(m,n))); 
            X=randn(m,n);
            imageflag=0;
            sparseflag=0;
        case 'rand-lr' % Random low-rank example
            m=16; 
            n=m; 
            r=floor(sqrt(min(m,n))); 
            X=rand(m,2*r)*rand(2*r,n);
            imageflag=0;
            sparseflag=0;
        case 'rand-hlr' % Random example with exact Hadamard Decomposition
            m=4096; 
            n=m;
            r1=12;
            r2=r1;
            r=[r1,r2];
            W1=randn(m,r1); 
            H1=randn(n,r1);
            W2=randn(m,r2);
            H2=randn(n,r2);
            X=(W1*H1').*(W2*H2');
            imageflag=0;
            sparseflag=0;
        case 'rand-sparse' % Random sparse example
            m=500;
            n=m;
            r=5;
            density=0.05;
            X=sprand(m,n,density);
            imageflag=0;
            sparseflag=1;
        case 'rand-sparse-hlr' % Random sparse Hadamard-Decomposable example
            m=25;
            n=m;
            r=5;
            density=0.3;
            beta=sqrt(density);
            [~,~,X1]=sparserandNMF(m,n,r,beta);
            [~,~,X2]=sparserandNMF(m,n,r,beta);
            X=sparse(X1).*sparse(X2);
            imageflag=0;
            sparseflag=1;
        case 'cameraman' % Cameraman image
            m=256;
            n=m;
            r=16;
            X=double(imread('cameraman.tif'));
            imageflag=1;
            sparseflag=0;
        case 'sparse_model'
            i=1;
            models={'NG20','classic','la12','ohscal','sports','la1',...
                'la2','reviews','k1b','hitech','tr41','tr45','tr11','tr23'};
            ranks=[20 4 6 10 7 6 6 5 6 6 10 10 9 6];
            string=['./datasets/',models{i},'.mat'];
            X=load(string).dtm;
            [m,n]=size(X);
            r=ranks(i);
        otherwise 
            error(['Example type not available. Try again with ',...
                '"rand-gen", "rand-lr", "rand-hlr", "rand-sparse" ',...
                '"rand-sparse-hlr", "cameraman" or "sparse_model".'])
    end

    %% 2) Select the parameters 
    % Choose optimization algorithm among 'Manopt','manBCD','projBCD' and
    % 'BCD'. Define the list methods={'method1,...,methodk}.
    methods={'Manopt','manBCD','projBCD','BCD'};

    % Choose initialization among 'SVD-based','FS','FSL' and 'FSR'. If you
    % want to select the best use 'all'. If you have given values for the
    % initialization, set them as opts.W1=W1, opts.H1=H1, opts.W2=W2 and
    % opts.H2=H2, by recalling that W1,W2 are m-by-r and H1,H2 are n-by-r.
    opts.init='FS';

    % Choose the time limit (in seconds) for each method in 'methods'
    opts.maxtime=200;

    % If your Matlab version is R2021b or newer, select opts.noloops=1,
    % (requires pagesvd built-in function), otherwise set opts.noloops=0.
    opts.noloops=1;


    %% 3) Main computation
    % Data storage
    n_methods=length(methods);
    W1=cell(n_methods,1);
    W2=cell(n_methods,1);
    H1=cell(n_methods,1);
    H2=cell(n_methods,1);
    relerr=zeros(n_methods+1,1);
    info=cell(n_methods+1,1);
    Xapprox=cell(n_methods+1,1);

    % Rank-r HDs
    for k=1:n_methods
        opts.method=methods{k};
        [W1{k},H1{k},W2{k},H2{k},info{k}]=HadDec(X,r,opts);
        relerr(k)=info{k}.err(end);
        Xapprox{k}=(W1{k}*H1{k}').*(W2{k}*H2{k}');
    end

    % TSVD of rank 2r
    rsvd=min([n,m,2*r]);
    normX=norm(X,'fro');
    Xsvd=X/normX;
    if sparseflag
        [Ursvd,Srsvd,Vrsvd]=svds(Xsvd,rsvd);
    else
        [U,S,V]=svd(Xsvd);
        Srsvd=S(1:rsvd,1:rsvd);
        Ursvd=U(:,1:rsvd);
        Vrsvd=V(:,1:rsvd);
    end
    Xapprox{end}=normX*Ursvd*Srsvd*Vrsvd';
    relerr(end)=norm(X-Xapprox{end},'fro')/normX;
   
    %% 4) Display the results
    fprintf('Final relative error(s): ')
    for k=1:n_methods
        fprintf('%s = %2.2f%%, ', methods{k}, relerr(k)*100);
    end
    fprintf('%i-%s = %2.2f%%\n', rsvd, 'TSVD', relerr(end)*100);

    % Plot relative error versus time
    close all
    lw=1.3;
    plotsettings={'m-','r-','b-','k-','g-'};
    figure
    for k=1:n_methods
        semilogy(info{k}.time,info{k}.err,plotsettings{k},'LineWidth',lw)
        hold on
    end
    f=relerr(end);
    maxt=0;
    for k=1:n_methods
        time_method=info{k}.time(end);
        if time_method>maxt
            maxt=time_method;
        end
    end
    semilogy([0,maxt],[f,f],'-.','Color',[0.75,0.5,0],'LineWidth',lw)
    hold on
    title('Relative error(s) vs time')
    legend([methods,{num2str(rsvd)+"-TSVD"}],'Location','best')

    %% For example='cameraman', plot the image(s) reconstructed 
    if imageflag
        figure
        switch n_methods
            case {1,2}
                a=2; b=2;
            case {3,4}
                a=2; b=3;
            case 5
                a=2; b=4;
        end
        subplot(a,b,1); 
        imshow(X/256); 
        title('Original image')
        for k=1:n_methods  
            subplot(a,b,k+1); 
            imshow(Xapprox{k}/256); 
            title([methods{k},': ',num2str(round(relerr(k)*100,2)),'%'])
            sgtitle('Relative errors of HD(s) and TSVD')
        end
        subplot(a,b,k+2); 
        imshow(Xapprox{end}/256); 
        title([num2str(rsvd),'-TSVD: ',num2str(round(relerr(end)*100,2)),'%'])
    end
