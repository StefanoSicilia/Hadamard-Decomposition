
    n=10;
    m=10;
    r=3;
    rng(1)
    Xtrue=randi(10,n,r);
    Ytrue=randi(10,m,r);
    Utrue=randi(10,n,r);
    Vtrue=randi(10,m,r);
    E=randi(10,n,m);
    epsilon=0;
    Atrue=(Xtrue*Ytrue').*(Utrue*Vtrue');
    A=Atrue/norm(Atrue,'fro')+epsilon*E/norm(E,'fro');

    opts=struct('r',r,'maxit',10,'init','FS',...
        'X',Xtrue,'Y',Ytrue,'U',Utrue,'V',Vtrue);
    [X_FS,Y_FS,U_FS,V_FS,err_FS]=HadBCD_alt(A,opts);
    fin_err_FS=err_FS(end);
    opts.init='given';
    [X_given,Y_given,U_given,V_given,err_given]=HadBCD_alt(A,opts);
    fin_err_given=err_given(end);

    close all
    lw=1.3;
    semilogy(err_FS,'r-','LineWidth',lw)
    hold on
    semilogy(err_given,'b-','LineWidth',lw)
    legend('FS','given','Location','best')
    