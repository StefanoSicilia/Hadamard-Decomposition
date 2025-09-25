function [W1,H1,W2,H2,err]=Had_manBCD(X,opts)
%% Had_manBCD:
% Computes a(n) (approximate) Hadamard decomposition X=WH', where 
% W=face_split(W1,W2) and H=face_split(H1,H2) have rank r=opts.rank. 
% It uses a 2 block coordinate descend algorithm for W and H and optimizes 
% each rank-r matrix on the manifold such that they admit a face-split 
% decomposition.

    [n,m]=size(X);
    r=opts.rank;
    maxit=opts.maxit;
    tol=opts.tol;
    extrapar=opts.momentum;
    beta=extrapar(1); betat=extrapar(2); gamma=extrapar(3); 
    gammat=extrapar(4); eta=extrapar(5);
    switch opts.init
        case 'Wertz'
            M=sqrt(abs(X));
            N=sign(X).*M;
            [U1,S1,V1]=svd(M);
            W1=U1(:,1:r)*sqrt(S1(1:r,1:r));
            H1=V1(:,1:r)*sqrt(S1(1:r,1:r));
            [U2,S2,V2]=svd(N);
            W2=U2(:,1:r)*sqrt(S2(1:r,1:r));
            H2=V2(:,1:r)*sqrt(S2(1:r,1:r));
        case 'FS'
            [W1,W2,H1,H2]=Had_init(X,r);
        case 'FS2'
            [W1,W2,H1,H2]=Had_init2(X,r);
        case 'FS3'
            [W1,W2,H1,H2]=Had_init3(X,r);
        case 'given'
            W1=opts.X;
            H1=opts.Y;
            W2=opts.U;
            H2=opts.V;
        otherwise
            error('Initialization not available.')
    end

    % Second initialization
    % W=face_split(W1,W2);
    % H=face_split(H1,H2);
    % G=(W*H'-X)*H;
    % for i=1:n
    %     Gi=reshape(G(i,:),r,r);
    %     [P,~,Q]=svd(Gi);
    %     W1(i,:)=P(1,:);
    %     W2(i,:)=Q(1,:);
    % end
    % G=(H*W'-A')*W;
    % for i=1:m
    %     Gi=reshape(G(i,:),r,r);
    %     [P,~,Q]=svd(Gi);
    %     H1(i,:)=P(1,:);
    %     H2(i,:)=Q(1,:);
    % end

    % normalization of the rows of X and Y
    for i=1:n
        normXi=norm(W1(i,:));
        W1(i,:)=W1(i,:)/normXi;
        W2(i,:)=W2(i,:)*normXi;
    end
    for i=1:m
        normYi=norm(H1(i,:));
        H1(i,:)=H1(i,:)/normYi;
        H2(i,:)=H2(i,:)*normYi;
    end
    W=face_split(W1,W2);
    H=face_split(H1,H2);
    normX=norm(X,'fro');
    
    Hblock=opts.Hblock;
    Wblock=opts.Wblock;
    err=zeros(maxit*(Hblock+Wblock)+1,1);
    err(1)=0.5*norm(X-W*H','fro')^2/normX;
    i=1;
    j=1;
    tau=opts.tau;

    if strcmp(opts.init,'FS3')  
        while j<=maxit && err(i)>tol
            % H-update
            for k=1:Hblock
                H1old=H1;
                H2old=H2;
                G=(H*W'-X')*W;
                %alpha=tau/norm(W'*W,'fro');
                alpha=tau/norm(W'*W,2);
                %[H1,H2]=Upd_manifold(H1,H2,pinv(W'*W)*G,alpha);
                [H1,H2]=Upd_manifold(H1,H2,G,alpha);
                H1=H1+beta*(H1-H1old);
                H2=H2+beta*(H2-H2old);
                H=face_split(H1,H2);
                i=i+1;
                err(i)=0.5*norm(X-W*H','fro')^2/normX;
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
                W1old=W1;
                W2old=W2;
                G=(W*H'-X)*H;
                %alpha=tau/norm(H'*H,'fro');
                alpha=tau/norm(H'*H,2);
                %[W1,W2]=Upd_manifold(W1,W2,pinv(H'*H)*G,alpha);
                [W1,W2]=Upd_manifold(W1,W2,G,alpha);
                W1=W1+beta*(W1-W1old);
                W2=W2+beta*(W2-W2old);
                W=face_split(W1,W2);
                i=i+1;
                err(i)=0.5*norm(X-W*H','fro')^2/normX;
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
    else
        while j<=maxit && err(i)>tol
            % W-update
            for k=1:Wblock
                W1old=W1;
                W2old=W2;
                G=(W*H'-X)*H;
                %alpha=tau/norm(H'*H,'fro');
                alpha=tau/norm(H'*H,2); 
                %[W1,W2]=Upd_manifold(W1,W2,pinv(H'*H)*G,alpha);
                [W1,W2]=Upd_manifold(W1,W2,G,alpha);
                W1=W1+beta*(W1-W1old);
                W2=W2+beta*(W2-W2old);
                W=face_split(W1,W2);
                i=i+1;
                err(i)=0.5*norm(X-W*H','fro')^2/normX;
                if err(i)<err(i-1)
                    beta=min(betat,gamma*beta);
                    betat=min(1,gammat*betat);
                else
                    beta=beta/eta;
                    betat=beta;
                end
            end

            % H-update
            for k=1:Hblock
                H1old=H1;
                H2old=H2;
                G=(H*W'-X')*W;
                %alpha=tau/norm(W'*W,'fro');
                alpha=tau/norm(W'*W,2);
                %[H1,H2]=Upd_manifold(H1,H2,pinv(W'*W)*G,alpha);
                [H1,H2]=Upd_manifold(H1,H2,G,alpha);
                H1=H1+beta*(H1-H1old);
                H2=H2+beta*(H2-H2old);
                H=face_split(H1,H2);
                i=i+1;
                err(i)=0.5*norm(X-W*H','fro')^2/normX;
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
    end

    err(i)=0.5*norm(X-W*H','fro')^2/normX;
    err=err(1:i);

end
