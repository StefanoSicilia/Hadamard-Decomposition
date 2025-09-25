function [W1,W2,H1,H2,err]=Had_eBCD(X,opts)
%% Had_eBCD:
% Code from the paper by Wertz et al. with extrapolation (momentum).
    
    init=opts.init;
    r=opts.rank;
    maxit=opts.maxit;
    extrapar=opts.momentum;
    beta=extrapar(1); betat=extrapar(2); gamma=extrapar(3); 
    gammat=extrapar(4); eta=extrapar(5);
    % [n,m]=size(A);
    % if min(n,m)<r^2
    %     error('The selected value of the rank is too high.')
    % end
    switch init
        case 'Wertz'
            M=sqrt(abs(X));
            N=sign(X).*M;
            [U1,S1,V1]=svd(M);
            W1=U1(:,1:r)*sqrt(S1(1:r,1:r));
            W2=V1(:,1:r)*sqrt(S1(1:r,1:r));
            [U2,S2,V2]=svd(N);
            H1=U2(:,1:r)*sqrt(S2(1:r,1:r));
            H2=V2(:,1:r)*sqrt(S2(1:r,1:r));
        case 'FS'
            [W1,H1,W2,H2]=Had_init(X,r);
        case 'FS2'
            [W1,H1,W2,H2]=Had_init2(X,r);
        case 'FS3'
            [W1,H1,W2,H2]=Had_init3(X,r);
        case 'given'
            W1=opts.X;
            W2=opts.Y;
            H1=opts.U;
            H2=opts.V;
        otherwise
            error('Initialization not available.')
    end

    % Orthogonalization of initialization
    % [Q1,R1]=qr(H1,0);
    % H1=Q1;
    % W1=W1*R1';
    % [Q2,R2]=qr(H2,0);
    % H2=Q2;
    % W2=W2*R2';

    err=zeros(maxit,1);
    normX=norm(X,'fro');
    err(1)=norm(X-(W1*W2').*(H1*H2'),'fro')/normX;
    i=1;
    while i<=maxit && err(i)>opts.tol
        i=i+1;
        H2old=H2; H1old=H1; W2old=W2; W1old=W1;
        H2=UpdFact(X,W1,W2,H1,H2);
        H2=H2+beta*(H2-H2old);
        H1=UpdFact(X',W2,W1,H2,H1);
        H1=H1+beta*(H1-H1old);
        W2=UpdFact(X,H1,H2,W1,W2);
        W2=W2+beta*(W2-W2old);
        W1=UpdFact(X',H2,H1,W2,W1);
        W1=W1+beta*(W1-W1old);
        err(i)=norm(X-(W1*W2').*(H1*H2'),'fro')/normX;
        if err(i)<err(i-1)
            beta=min(betat,gamma*beta);
            betat=min(1,gammat*betat);
        else
            beta=beta/eta;
            betat=beta;
        end
    end
    err(i)=norm(X-(W1*W2').*(H1*H2'),'fro')/normX;
    err=err(1:i);
end