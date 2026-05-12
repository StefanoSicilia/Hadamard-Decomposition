%% Script to test HadDec on simple examples for rank-r HD
% To run the test, follow the steps and choose among the possible options.
% The script takes 20 seconds to run (in general opts.maxtime*n_methods).

    %% 1) Choose an example: 'rand-gen','rand-lr','rand-hlr' or 'cameraman'
    example='cameraman';
    
    rng(1)
    switch example
        case 'rand-gen' % Random generic example
            m=16; 
            n=m; 
            r=floor(sqrt(min(m,n))); 
            X=rand(m,n);
            imageflag=0;
        case 'rand-lr' % Random low-rank example
            m=16; 
            n=m; 
            r=floor(sqrt(min(m,n))); 
            X=rand(m,2*r)*rand(2*r,n);
            imageflag=0;
        case 'rand-hlr' % Random example with exact Hadamard Decomposition
            m=16; 
            n=m; 
            r=floor(sqrt(min(m,n))); 
            W1=randi(m,r1); 
            H1=randi(n,r1);
            W2=randi(m,r2);
            H2=randi(n,r2);
            X=(W1*H1').*(W2*H2');
            imageflag=0;
        case 'cameraman' % Cameraman image
            m=256;
            n=m;
            r=16;
            X=double(imread('cameraman.tif'));
            imageflag=1;
        otherwise 
            error(['Example type not available. Try again with ',...
                '"rand-gen", "rand-lr", "rand-hlr" or "cameraman".'])
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
    opts.maxtime=5;

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
    [U,S,V]=svd(Xsvd);
    Srsvd=S(1:rsvd,1:rsvd);
    Ursvd=U(:,1:rsvd);
    Vrsvd=V(:,1:rsvd);
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
    plotsettings={'m-','r-','b-','k-'};
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
