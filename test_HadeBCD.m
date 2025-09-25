
    n=16;
    m=n;
    r=sqrt(n);
    nmax=3;
    nsample=10;
    maxit=500;

    rng(11)
    Xtrue=randi(nmax,n,r);
    Ytrue=randi(nmax,m,r);
    Utrue=randi(nmax,n,r);
    Vtrue=randi(nmax,m,r);
    Atrue=(Xtrue*Ytrue').*(Utrue*Vtrue');
    %Atrue=randi(nmax,n,2*r)*randi(nmax,2*r,m);
    %Atrue=randi(nmax,n,r); Atrue=Atrue*Atrue';
    %Atrue=ones(n,m)-eye(n,m);
    %Atrue=randi(n,m);
    %Atrue(1,2)=0; 
    %Atrue(2,3)=-1; Atrue(4,1)=10; Atrue(2,1)=0; Atrue(1,3)=0;

    opts=struct('rank',r,'maxit',maxit,'init','FS',...
        'momentum',[0.75,1,1.05,1.01,1.5]);
    [X_FS,Y_FS,U_FS,V_FS,err_FS]=Had_eBCD(Atrue,opts);
    fin_err_FS=err_FS(end);
    opts.init='Wertz';
    [X_Wertz,Y_Wertz,U_Wertz,V_Wertz,err_Wertz]=Had_eBCD(Atrue,opts);
    fin_err_Wertz=err_Wertz(end);
    opts.init='FS2';
    tic;
    [X_FS2,Y_FS2,U_FS2,V_FS2,err_FS2]=Had_eBCD(Atrue,opts);
    t_Sam=toc;
    fin_err_FS2=err_FS2(end);
    opts.init='FS3';
    [X_FS3,Y_FS3,U_FS3,V_FS3,err_FS3]=Had_eBCD(Atrue,opts);
    fin_err_FS3=err_FS3(end);

    X_rand=zeros(n,r,nsample);
    Y_rand=zeros(m,r,nsample);
    U_rand=zeros(n,r,nsample);
    V_rand=zeros(m,r,nsample);
    err_rand=zeros(maxit,nsample);
    fin_err_rand=zeros(nsample,1);
    opts.init='given';
    for j=1:nsample
        rng(j)
        opts.X=randi(nmax,n,r); opts.U=randi(nmax,n,r);
        opts.Y=randi(nmax,m,r); opts.V=randi(nmax,m,r);
        [X_rand(:,:,j),Y_rand(:,:,j),U_rand(:,:,j),V_rand(:,:,j),...
            err_rand(:,j)]=Had_eBCD(Atrue,opts);
        fin_err_rand(j)=err_rand(end,j);
    end

    close all
    lw=1.3;
    semilogy(err_FS,'r-','LineWidth',lw)
    hold on
    semilogy(err_Wertz,'b-','LineWidth',lw)
    hold on
    semilogy(err_FS2,'g-','LineWidth',lw)
    hold on
    semilogy(err_FS3,'m-','LineWidth',lw)
    legend('FS','Wertz','FS2','FS3','Location','best')
    %legend('Wertz')

    figure
    semilogy(fin_err_rand,'c-','LineWidth',lw)
    title('Random samples error')
