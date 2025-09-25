function [X,U,Y,V,err]=Had_riemopt(A,opts)
%% Had_manopt:
% Computes a(n) (approximate) Hadamard decomposition A=WH', where 
% W=face_split(X,U) and H=face_split(Y,V) have rank r=opts.rank. 
% It uses a 2 block coordinate descend algorithm for W and H and optimizes 
% each rank-r matrix through Manopt.

    [n,m]=size(A);
    r=opts.rank;
    maxit=opts.maxit;
    h=opts.h;
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
    % normalization of the rows of X and Y
    for i=1:n
        normXi=norm(X(i,:));
        X(i,:)=X(i,:)/normXi;
        U(i,:)=U(i,:)*normXi;
    end
    for i=1:m
        normYi=norm(Y(i,:));
        Y(i,:)=Y(i,:)/normYi;
        V(i,:)=V(i,:)*normYi;
    end
    W=face_split(X,U);
    H=face_split(Y,V);
    err=zeros(2*maxit,1);

    alpha=opts.alpha;
    
    if strcmp(opts.init,'FS3')  
        for j=1:maxit
            % H-update
            G=2*(H*W'-A')*W;
            % InvHess=.5*(pinv(W'*W));
            % G=InvHess*G;
            for i=1:m
                Gi=reshape(G(i,:),r,r);
                yi=Y(i,:)';
                vi=V(i,:)';
                Y(i,:)=Y(i,:)+h*(-Gi'*vi+(vi'*Gi*yi)*yi)';
                Y(i,:)=Y(i,:)/norm(Y(i,:));
                V(i,:)=V(i,:)+h*(-Gi*yi)';
            end
            H=face_split(Y,V);
            err(2*j-1)=norm(A-W*H','fro')^2;
            h=alpha*err(2*j-1)/norm(G,'fro');
            
            % W-update
            G=2*(W*H'-A)*H;
            % InvHess=.5*(pinv(H'*H));
            % G=InvHess*G;
            for i=1:n
                Gi=reshape(G(i,:),r,r);
                xi=X(i,:)';
                ui=U(i,:)';
                X(i,:)=X(i,:)+h*(-Gi'*ui+(ui'*Gi*xi)*xi)';
                X(i,:)=X(i,:)/norm(X(i,:));
                U(i,:)=U(i,:)+h*(-Gi*xi)';
            end
            W=face_split(X,U);
            err(2*j)=norm(A-W*H','fro')^2;
            h=alpha*err(2*j)/norm(G,'fro');
        end
    else
        for j=1:maxit
            % W-update
            G=2*(W*H'-A)*H;
            % InvHess=.5*(pinv(H'*H));
            % G=InvHess*G;
            for i=1:n
                Gi=reshape(G(i,:),r,r);
                xi=X(i,:)';
                ui=U(i,:)';
                X(i,:)=X(i,:)+h*(-Gi'*ui+(ui'*Gi*xi)*xi)';
                X(i,:)=X(i,:)/norm(X(i,:));
                U(i,:)=U(i,:)+h*(-Gi*xi)';
            end
            W=face_split(X,U);
            err(2*j-1)=norm(A-W*H','fro')^2;
            h=alpha*err(2*j-1)/norm(G,'fro');

            % H-update
            G=2*(H*W'-A')*W;
            % InvHess=.5*(pinv(W'*W));
            % G=InvHess*G;
            for i=1:m
                Gi=reshape(G(i,:),r,r);
                yi=Y(i,:)';
                vi=V(i,:)';
                Y(i,:)=Y(i,:)+h*(-Gi'*vi+(vi'*Gi*yi)*yi)';
                Y(i,:)=Y(i,:)/norm(Y(i,:));
                V(i,:)=V(i,:)+h*(-Gi*yi)';
            end
            H=face_split(Y,V);
            err(2*j)=norm(A-W*H','fro')^2;
            h=alpha*err(2*j)/norm(G,'fro');
        end
    end

end
