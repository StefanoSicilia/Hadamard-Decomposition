%% HadDec on dog colored image
% It takes a bit less than 1 hour.
% The computational time is approximately given by 
% 3*n_methods*opts.maxtime

    %% Methods, parameters and structures
    m=400;
    n=600;
    r=20;
    X=double(imread('./datasets/rgb_dog.jpg'));
    X(m,:,:)=X(m-1,:,:);

    %% Methods parameters
    maxit=1e6;
    tol=1e-16;
    Iter_W=3;
    Iter_H=3;
    rng(1)
    opts=struct('maxit',maxit,'init','all','tau',0.95,'theta',1e-4,...
        'Iter_W',Iter_W,'Iter_H',Iter_H,'tol',tol,'sparsity',1,'noloop',1); 
    opts.momentum=[0.75,1,1.05,1.01,1.5,0.6];
    opts.maxtime=240; 
    r2=min([m,n,2*r]);
    methods={'Manopt','manBCD','projBCD','BCD','TSVD'};
    n_methods=length(methods)-1;
    W1=cell(3,n_methods);
    W2=cell(3,n_methods);
    H1=cell(3,n_methods);
    H2=cell(3,n_methods);
    times=zeros(3,n_methods+1);
    relerr=zeros(3,n_methods+1);
    info=cell(3,n_methods+1);
    normX=zeros(3,1);
    Xsvd=zeros(size(X));

    U=cell(3,1);
    S=cell(3,1);
    V=cell(3,1);
    Ur2=cell(3,1);
    Sr2=cell(3,1);
    Vr2=cell(3,1);

    %% Rank-r HDs and TSVD of rank 2r of each slice of the image
    for k=1:3
        
        for j=1:n_methods
            opts.method=methods{j};
            tic;
            [W1{k,j},H1{k,j},W2{k,j},H2{k,j},info{k,j}]=HadDec(X(:,:,k),r,opts);
            times(k,j)=toc;
            relerr(k,j)=info{k,j}.err(end);
        end
        
        normX(k)=norm(X(:,:,k),'fro');
        Xsvd(:,:,k)=X(:,:,k)/normX(k);
        [U{k},S{k},V{k}]=svd(Xsvd(:,:,k));
        Sr2{k}=S{k}(1:r2,1:r2);
        Ur2{k}=U{k}(:,1:r2);
        Vr2{k}=V{k}(:,1:r2);
        relerr(k,end)=norm(Xsvd(:,:,k)-Ur2{k}*Sr2{k}*Vr2{k}','fro');

    end


    %% Display the relative errors
    close all
    lw=1.3;   
    plotsettings={'m-','r-','b-','k-'};
    for k=1:3
        figure
        for j=1:n_methods
            semilogy(info{k,j}.time,info{k,j}.err,plotsettings{j},'LineWidth',lw)
            hold on
        end
        xlabel('Time (s)')
        ylabel('Objective function')
        title(['Slice ',num2str(k)])
        f=relerr(end);
        maxt=0;
        for j=1:n_methods
            time_method=info{k,j}.time(end);
            if time_method>maxt
                maxt=time_method;
            end
        end
        semilogy([0,maxt],[f,f],'-.','Color',[0.75,0.5,0],'LineWidth',lw)
        hold on
        legend(methods,'Location','best')
    end

    %% Plot the image compressions
    sq=255;
    A=cell(n_methods+1);
    for j=1:n_methods
        figure
        for k=1:3
            A{j}(:,:,k)=((W1{k,j}*H1{k,j}').*(W2{k,j}*H2{k,j}'))/normX(k);
        end
        imshow(sq*A{j})
        title(methods{j})
    end

    figure
    for k=1:3
        A{end}(:,:,k)=(Ur2{k}*Sr2{k}*Vr2{k}');
    end
    imshow(sq*A{end})
    title(methods{end})
    figure
    imshow(sq*Xsvd)
    title('original')


    %% Initializations chosen and number of iterations
    Inits=cell(3,n_methods);
    Iters=zeros(3,n_methods);
    for k=1:3
        for j=1:n_methods
            Inits{k,j}=info{k,j}.init;
            Iters(k,j)=length(info{k,j}.err);
        end
    end