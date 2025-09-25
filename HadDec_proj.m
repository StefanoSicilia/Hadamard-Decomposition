function [W,H,err]=HadDec_proj(A,opts)
%% HadDec_proj:
% Computes a(n) (approximate) Hadamard decomposition A=WH', where 
% W=face_split(X,U) and H=face_split(Y,V) have rank r=opts.rank. 
% It uses a 2 block coordinate descend algorithm for W and H and then
% projects onto the manifold to which W abd H belong.

    r=opts.rank;
    maxit=opts.maxit;
    switch opts.init
        case 'Wertz'
            M=sqrt(abs(A));
            N=sign(A).*M;
            [U1,S1,V1]=svd(M);
            X=U1(:,1:r)*sqrt(S1(1:r,1:r));
            Y=V1(:,1:r)*sqrt(S1(1:r,1:r));
            [U2,S2,V2]=svd(N);
            U=U2(:,1:r)*sqrt(S2(1:r,1:r));
            V=V2(:,1:r)*sqrt(S2(1:r,1:r));
            W=face_split(X,U);
            H=face_split(Y,V);
        case 'FS'
            [X,U,Y,V]=Had_init(A,r);
            W=face_split(X,U);
            H=face_split(Y,V);
        case 'FS2'
            [X,U,Y,V]=Had_init2(A,r);
            W=face_split(X,U);
            H=face_split(Y,V);
        case 'FS3'
            [X,U,Y,V]=Had_init3(A,r);
            W=face_split(X,U);
            H=face_split(Y,V);
        case 'given'
            W=opts.W;
            H=opts.H;
        otherwise
            error('Initialization not available.')
    end
   
    
    err=zeros(2*maxit,1);
    err(1)=norm(A-W*H','fro')^2;
    i=2;

    for j=1:maxit
        % W update
        W=ProjRowRK1(A*pinv(H)');
        err(i)=norm(A-W*H','fro')^2;
        i=i+1;

        % H update
        H=ProjRowRK1(A'*pinv(W));
        err(i)=norm(A-W*H','fro')^2;
        i=i+1;
    end

    err(i)=norm(A-W*H','fro')^2;
    err=err(1:i);

end
