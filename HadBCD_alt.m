function [X,U,Y,V,err]=HadBCD_alt(A,opts)

    init=opts.init;
    r=opts.r;
    [n,m]=size(A);
    if min(n,m)<r^2
        error('The selected value of the rank is too high.')
    end
    switch init
        case 'FS'
            [X,U,Y,V]=Had_init(A,r);
        case 'given'
            Y=opts.Y;
            V=opts.V;
            U=opts.U;
        otherwise
            error('Initialization not available.')
    end
    err=zeros(opts.maxit,1);
    normA=norm(A,'fro');
    for i=1:opts.maxit
        X=UpdFact_alt(A,Y,V,U);
        Y=UpdFact_alt(A',X,U,V);
        U=UpdFact_alt2(A,Y,V,X);
        V=UpdFact_alt2(A',X,U,Y);
        err(i)=norm(A-(X*Y').*(U*V'),'fro')/normA;
    end


end