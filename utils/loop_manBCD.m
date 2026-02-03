function [W1,H1,W2,H2,info]=loop_manBCD(X,W1,H1,W2,H2,opts)
%% loop_manBCD: Main cycle for HadDec - manBCD case
% Performs Riemannian gradient descent for Hadamard decomposition. See
% HadDec_manBCD for more details.

    % Initialization and parameters
    tstart=tic;
    [m,n]=size(X);
    maxit=opts.maxit;
    tol=opts.tol;
    extrapar=opts.momentum;
    beta=extrapar(1); betat=extrapar(2); gamma=extrapar(3); 
    gammat=extrapar(4); eta=extrapar(5);
    tau=opts.tau;

    % Normalization of the rows of W1 and H1
    % for i=1:m
    %     normXi=norm(W1(i,:));
    %     W1(i,:)=W1(i,:)/normXi;
    %     W2(i,:)=W2(i,:)*normXi;
    % end
    % for i=1:n
    %     normYi=norm(H1(i,:));
    %     H1(i,:)=H1(i,:)/normYi;
    %     H2(i,:)=H2(i,:)*normYi;
    % end
    W=face_split(W1,W2);
    H=face_split(H1,H2);

    Hblock=opts.Hblock;
    Wblock=opts.Wblock;
    err=zeros(maxit,1);
    err(1)=0.5*norm(X-W*H','fro')^2;
    i=1;
    j=1;
    time=zeros(maxit,1);
    init_time=toc(tstart);
    time(1)=init_time;
    maxtime=opts.maxtime-init_time;
    tstart=tic;
    

    % Main cycle
    while j<maxit && err(i)>tol && time(i)<maxtime
        % W-update
        for k=1:Wblock
            W1old=W1;
            W2old=W2;
            G=(W*H'-X)*H;
            alpha=tau/norm(H'*H,2); 
            [W1,W2]=Upd_manifold(W1,W2,G,alpha);
            W1=W1+beta*(W1-W1old);
            W2=W2+beta*(W2-W2old);
            W=face_split(W1,W2);
            i=i+1;
            err(i)=0.5*norm(X-W*H','fro')^2;
            time(i)=toc(tstart);
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
            alpha=tau/norm(W'*W,2);
            [H1,H2]=Upd_manifold(H1,H2,G,alpha);
            H1=H1+beta*(H1-H1old);
            H2=H2+beta*(H2-H2old);
            H=face_split(H1,H2);
            i=i+1;
            err(i)=0.5*norm(X-W*H','fro')^2;
            time(i)=toc(tstart);
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

    % Output
    err(i)=0.5*norm(X-W*H','fro')^2;
    info=struct('err',sqrt(2*err(1:i)),'time',time(1:i));


end

