function [X,U,Y,V,err]=HadDec_manifold(A,opts)
%% Had_manopt:
% Computes a(n) (approximate) Hadamard decomposition A=WH', where 
% W=face_split(X,U) and H=face_split(Y,V) have rank r=opts.rank. 
% It uses a 2 block coordinate descend algorithm for W and H and optimizes 
% each rank-r matrix through Manopt.

    [n,m]=size(A);
    r=opts.rank;
    maxit=opts.maxit;
    theta=opts.theta;
    safestop=opts.safestop;
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
    

    c=opts.c;
    beta=opts.beta;
    Hblock=opts.Hblock;
    Wblock=opts.Wblock;
    err=zeros(maxit*(Hblock+Wblock)+1,1);
    f=norm(A-W*H','fro')^2;
    err(i)=f;
    i=2;
    tau=opts.tau;

    if strcmp(opts.init,'FS3')  
        for j=1:maxit
            % H-update
            Hold=H;
            for k=1:Hblock
                G=2*(H*W'-A')*W;
                fh=f;
                normG=norm(G,'fro');
                %alpha=tau*f/normG;
                alpha=0.5/norm(W'*W,'fro');
                cont=1;
                while fh>=f-alpha*c*normG^2 && cont<safestop
                    alpha=alpha/theta;
                    [Y,V]=Upd_manifold(Y,V,G,alpha);
                    H=face_split(Y,V);
                    fh=norm(A-W*H','fro')^2;
                    cont=cont+1;
                    f=fh;
                end
                err(i)=f;
                i=i+1;
                f=fh;
            end
            H=H+beta*(H-Hold);
            
            % W-update
            Wold=W;
            for k=1:Wblock
                G=2*(W*H'-A)*H;
                fh=f;
                normG=norm(G,'fro');
                %alpha=tau*f/normG;
                alpha=0.5/norm(H'*H,'fro');
                cont=1;
                while fh>=f-alpha*c*normG^2 && cont<safestop
                    alpha=alpha/theta;
                    [X,U]=Upd_manifold(X,U,G,alpha);
                    W=face_split(X,U);
                    fh=norm(A-W*H','fro')^2;
                    cont=cont+1;
                    f=fh;
                end
                err(i)=f;
                i=i+1;
                f=fh;
            end
            W=W+beta*(W-Wold);
        end
    else
        for j=1:maxit
            % W-update
            Wold=W;
            for k=1:Wblock
                G=2*(W*H'-A)*H;
                fh=f;
                normG=norm(G,'fro');
                %alpha=tau*f/normG;
                alpha=0.5/norm(H'*H,'fro');
                cont=1;
                while fh>=f-alpha*c*normG^2 && cont<safestop
                    alpha=alpha/theta;
                    [X,U]=Upd_manifold(X,U,G,alpha);
                    W=face_split(X,U);
                    fh=norm(A-W*H','fro')^2;
                    cont=cont+1;
                    f=fh;
                end
                err(i)=f;
                i=i+1;
                f=fh;
            end
            W=W+beta*(W-Wold);

            % H-update
            Hold=H;
            for k=1:Hblock
                G=2*(H*W'-A')*W;
                fh=f;
                normG=norm(G,'fro');
                %alpha=tau*f/normG;
                alpha=0.5/norm(W'*W,'fro');
                cont=1;
                while fh>=f-alpha*c*normG^2 && cont<safestop
                    alpha=alpha/theta;
                    [Y,V]=Upd_manifold(Y,V,G,alpha);
                    H=face_split(Y,V);
                    fh=norm(A-W*H','fro')^2;
                    cont=cont+1;
                    f=fh;
                end
                err(i)=f;
                i=i+1;
                f=fh;
            end
            H=H+beta*(H-Hold);
        end
    end

    G1=2*(W*H'-A)*H;
    G2=2*(H*W'-A')*W;
    norm(G1,'fro')
    norm(G2,'fro')
    err(i)=norm(A-W*H','fro')^2;
    err=err(1:i);

end
