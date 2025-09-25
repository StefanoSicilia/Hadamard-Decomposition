function [X,Y,U,V,err]=Had_BCD(A,opts)
%% Had_BCD:
% Code from the paper by Wertz et al.
    
    init=opts.init;
    r=opts.rank;
    % [n,m]=size(A);
    % if min(n,m)<r^2
    %     error('The selected value of the rank is too high.')
    % end
    switch init
        case 'Wertz'
            M=sqrt(abs(A));
            N=sign(A).*M;
            [U1,S1,V1]=svd(M);
            X=U1(:,1:r)*sqrt(S1(1:r,1:r));
            Y=V1(:,1:r)*sqrt(S1(1:r,1:r));
            [U2,S2,V2]=svd(N);
            U=U2(:,1:r)*sqrt(S2(1:r,1:r));
            V=V2(:,1:r)*sqrt(S2(1:r,1:r));
        case 'FS'
            [X,U,Y,V]=Had_init(A,r);
        case 'FS2'
            [X,U,Y,V]=Had_init2(A,r);
        case 'FS3'
            [X,U,Y,V]=Had_init3(A,r);
        case 'given'
            X=opts.X;
            Y=opts.Y;
            U=opts.U;
            V=opts.V;
        otherwise
            error('Initialization not available.')
    end
    err=zeros(opts.maxit,1);
    normA=norm(A,'fro');
    for i=1:opts.maxit
        V=UpdFact(A,X,Y,U,V);
        U=UpdFact(A',Y,X,V,U);
        Y=UpdFact(A,U,V,X,Y);
        X=UpdFact(A',V,U,Y,X);
        err(i)=norm(A-(X*Y').*(U*V'),'fro')/normA;
    end
end