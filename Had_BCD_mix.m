function [X,Y,U,V,err]=Had_BCD_mix(A,opts)
%% Had_manBCD:
% Computes a(n) (approximate) Hadamard decomposition A=WH', where 
% W=face_split(X,U) and H=face_split(Y,V) have rank r=opts.rank. 
% It uses a 2 block coordinate descend algorithm for W and H and optimizes 
% each rank-r matrix on the manifold such that they admit a face-split 
% decomposition.

    [n,m]=size(A);
    r=opts.rank;
    maxit=opts.maxit;
    tol=opts.tol;
    extrapar=opts.momentum;
    beta=extrapar(1); betat=extrapar(2); gamma=extrapar(3); 
    gammat=extrapar(4); eta=extrapar(5);
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

    % Second initialization
    % W=face_split(X,U);
    % H=face_split(Y,V);
    % G=(W*H'-A)*H;
    % for i=1:n
    %     Gi=reshape(G(i,:),r,r);
    %     [P,~,Q]=svd(Gi);
    %     X(i,:)=P(1,:);
    %     U(i,:)=Q(1,:);
    % end
    % G=(H*W'-A')*W;
    % for i=1:m
    %     Gi=reshape(G(i,:),r,r);
    %     [P,~,Q]=svd(Gi);
    %     Y(i,:)=P(1,:);
    %     V(i,:)=Q(1,:);
    % end

    % normalization of the rows of X
    for i=1:n
        normXi=norm(X(i,:));
        X(i,:)=X(i,:)/normXi;
        U(i,:)=U(i,:)*normXi;
    end
    W=face_split(X,U);
    H=face_split(Y,V);
    normA=norm(A,'fro');
    
    Hblock=opts.Hblock;
    Wblock=opts.Wblock;
    err=zeros(maxit*(Hblock+Wblock)+1,1);
    err(1)=0.5*norm(A-W*H','fro')^2/normA;
    i=1;
    j=1;
    tau=opts.tau;

    while j<=maxit && err(i)>tol
        % H-update
        for k=1:Hblock
            Vold=V; Yold=Y; 
            V=UpdFact(A,X,Y,U,V);
            V=V+beta*(V-Vold);
            Y=UpdFact(A,U,V,X,Y);
            Y=Y+beta*(Y-Yold);
            H=face_split(Y,V);
            i=i+1;
            err(i)=0.5*norm(A-W*H','fro')^2/normA;
            if err(i)<err(i-1)
                beta=min(betat,gamma*beta);
                betat=min(1,gammat*betat);
            else
                beta=beta/eta;
                betat=beta;
            end
        end
        
        % W-update
        for k=1:Wblock
            Xold=X;
            Uold=U;
            G=(W*H'-A)*H;
            %alpha=tau/norm(H'*H,'fro');
            alpha=tau/norm(H'*H,2);
            %[X,U]=Upd_manifold(X,U,pinv(H'*H)*G,alpha);
            [X,U]=Upd_manifold(X,U,G,alpha);
            X=X+beta*(X-Xold);
            U=U+beta*(U-Uold);
            W=face_split(X,U);
            i=i+1;
            err(i)=0.5*norm(A-W*H','fro')^2/normA;
            if err(i)<err(i-1)
                beta=min(betat,gamma*beta);
                betat=min(1,gammat*betat);
            else
                beta=beta/eta;
                betat=beta;
            end
        end
        j=j+1;
    end

    err(i)=0.5*norm(A-W*H','fro')^2/normA;
    err=err(1:i);

end
