%% HadDec on dog colored image
% It takes about 30 minutes

    %% Methods, parameters and structures
    m=399;
    n=600;
    r=20;
    X=double(imread('./datasets/rgb_dog.jpg'));
    X(m+1,:,:)=ones(1,n,3);
    imageflag=1;

    maxit=1e6;
    tol=1e-15;
    Hblock=1;
    Wblock=1;
    rng(1)
    opts=struct('rank',r,'maxit',maxit,'init','all','tau',1.5,...
        'Hblock',Hblock,'Wblock',Wblock,'tol',0.5*tol^2);   
    opts.momentum=[0.75,1,1.05,1.01,1.5];
    %opts.momentum=[0,0,0,0,1];
    opts.maxtime=240; 
    r2=min([m,n,2*r]);

    % Note: total time required by HadDec is 9*opts.maxtime 

    W1=cell(3,3);
    W2=cell(3,3);
    H1=cell(3,3);
    H2=cell(3,3);
    info=cell(3,3);
    times=zeros(3);
    fin_err=zeros(3,4);
    normX=zeros(3,1);
    Xsvd=zeros(size(X));

    U=cell(3,1);
    S=cell(3,1);
    V=cell(3,1);
    Ur2=cell(3,1);
    Sr2=cell(3,1);
    Vr2=cell(3,1);

    opts_1=opts; opts_1.method='manBCD';
    opts_2=opts; opts_2.method='BCD';
    opts_3=opts; opts_3.method='Manopt';

    for k=1:3

        % manBCD
        j=1;
        [W1{k,j},H1{k,j},W2{k,j},H2{k,j},info{k,j}]=HadDec(X(:,:,k),r,opts_1);
        fin_err(k,j)=info{k,j}.err(end);

        % BCD
        j=2;
        [W1{k,j},H1{k,j},W2{k,j},H2{k,j},info{k,j}]=HadDec(X(:,:,k),r,opts_2);
        fin_err(k,j)=info{k,j}.err(end);

        % Manopt
        j=3;
        [W1{k,j},H1{k,j},W2{k,j},H2{k,j},info{k,j}]=HadDec(X(:,:,k),r,opts_3);
        fin_err(k,j)=info{k,j}.err(end);
        
        % SVD
        j=4;
        normX(k)=norm(X(:,:,k),'fro');
        Xsvd(:,:,k)=X(:,:,k)/normX(k);
        [U{k},S{k},V{k}]=svd(Xsvd(:,:,k));
        Sr2{k}=S{k}(1:r2,1:r2);
        Ur2{k}=U{k}(:,1:r2);
        Vr2{k}=V{k}(:,1:r2);
        fin_err(k,j)=norm(Xsvd(:,:,k)-Ur2{k}*Sr2{k}*Vr2{k}','fro');

    end

    %% Display the relative errors
    close all
    countfig=0;
    lw=1.3;   
    legendlabel={'manBCD','BCD','Manopt','2r-SVD'};
    colorlabel={'r','b','g',[0.75,0.5,0]};
    for k=1:3
        countfig=countfig+1;
        figure(countfig)
        for j=1:3
            semilogy(info{k,j}.time,info{k,j}.err,'-','LineWidth',lw,...
                'Color',colorlabel{j})
            hold on
        end
        f=fin_err(k,4);
        Mtime=max([info{k,1}.time(end),info{k,2}.time(end),info{k,3}.time(end)]);
        semilogy([0,Mtime],[f,f],'-.','Color',[0.75,0.5,0],'LineWidth',lw)
        legend(legendlabel,'Location','best')
    end

    %% Plot the image compressions
    sq=255;
    A=cell(4);
    title_label={'manBCD','BCD','Manopt','2r-SVD'};
    for j=1:3
        countfig=countfig+1;
        figure(countfig)
        for k=1:3
            A{j}(:,:,k)=((W1{k,j}*H1{k,j}').*(W2{k,j}*H2{k,j}'))/normX(k);
        end
        imshow(sq*A{j})
        title(title_label{j})
    end
    j=4;
    countfig=countfig+1;
    figure(countfig)
    for k=1:3
        A{j}(:,:,k)=(Ur2{k}*Sr2{k}*Vr2{k}');
    end
    imshow(sq*A{j})
    title(title_label{j})
    countfig=countfig+1;
    figure(countfig)
    imshow(sq*Xsvd)
    title('original')

    %% Slices relative errors
    for i=1:3
        figure(i)
        xlabel('Time (s)')
        ylabel('Objective function')
        title(['Slice ',num2str(i)])
    end

    %% Initializations chosen and number of iterations
    Inits=cell(3,3);
    Iters=zeros(3,3);
    for k=1:3
        for j=1:3
            Inits{k,j}=info{k,j}.init;
            Iters(k,j)=length(info{k,j}.err);
        end
    end