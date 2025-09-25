function [X,U,Y,V,err]=Had_manopt(A,opts)
%% Had_manopt:
% Computes a(n) (approximate) Hadamard decomposition A=WH', where 
% W=face_split(X,U) and H=face_split(Y,V) have rank r=opts.rank. 
% It uses a 2 block coordinate descend algorithm for W and H and optimizes 
% each rank-r matrix through Manopt.

    [n,m]=size(A);
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
    Mr=fixedrankembeddedfactory(r,r,1);
    manifoldW.X=powermanifold(Mr,n);
    manifoldW.U=powermanifold(Mr,n);
    manifoldW=productmanifold(manifoldW);
    manifoldH.Y=powermanifold(Mr,m);
    manifoldH.V=powermanifold(Mr,m);
    manifoldH=productmanifold(manifoldH);
    problemW.M=manifoldW;
    problemH.M=manifoldH;
    RW={cell(n,1),cell(n,1)};
    RH={cell(m,1),cell(m,1)};
    for i=1:n
        RW{1}{i}=X(i,:);
        RW{2}{i}=U(i,:);
    end
    for i=1:m
        RH{1}{i}=Y(i,:);
        RH{2}{i}=V(i,:);
    end
    W=face_split(X,U);
    H=face_split(Y,V);
    err=zeros(2*maxit,1);
    
    for j=1:maxit

        % W-update
        problemW.cost=@(R) resW(R,A,H);
        problemW.egrad=@(R) resgradW(R,A,H);
        problemW.ehess=@(R) 2*(H'*H);
        [RW, rescost, ~] = trustregions(problemW,RW);
        err(2*j-1)=rescost;
        for i=1:n
            X(i,:)=RW{1}{i};
            U(i,:)=RW{2}{i};
        end
        W=face_split(X,U);

        % H-update
        problemH.cost=@(R) resW(R,A',W);
        problemH.egrad=@(R) resgradW(R,A',W);
        problemH.ehess=@(R) 2*(W'*W);
        [RH,rescost, ~] = trustregions(problemW,RH);
        err(2*j)=rescost;
        for i=1:m
            Y(i,:)=RH{1}{i};
            V(i,:)=RH{2}{i};
        end
        H=face_split(Y,V);
    end

    % output reconstruction
    for i=1:n
        X(i,:)=RW{1}{i};
        U(i,:)=RW{2}{i};
    end
    for i=1:m
        Y(i,:)=RH{1}{i};
        V(i,:)=RH{2}{i};
    end

end
