function [W1,H1,W2,H2,info]=loop_manBCDres(X,W1,H1,W2,H2,opts)
%% loop_manBCDres: Main cycle for HadDec - manBCD case with restart
% Performs Riemannian gradient block coordinate descent for Hadamard 
% decomposition of X, with rank provided by the starting approximation 
% X~=(W1*H1').*(W2*H2'). The iterations computed are directly updated onto
% the manifold Bmr. See HadDec for more details.

    % Initialization and parameters
    t_start=tic;
    maxit=opts.maxit;
    tol=opts.tol;
    extrapar=opts.momentum;
    beta=extrapar(1); betat=extrapar(2); gamma=extrapar(3); 
    gammat=extrapar(4); eta=extrapar(5); betaold=beta;
    tau=opts.tau;
    Iter_W=opts.Iter_W;
    Iter_H=opts.Iter_H;
    j=1;
    time=zeros(maxit,1);
    init_time=toc(t_start);
    err=zeros(maxit,1);
    err(j)=norm(X-face_split(W1,W2)*face_split(H1,H2)','fro');  
    time(j)=init_time;
    maxtime=opts.maxtime-init_time;
    W1y=W1;
    W2y=W2;
    H1y=H1;
    H2y=H2;

    % Main cycle
    while j<maxit && err(j)>tol && time(j)<maxtime

        t_iter=tic;
        [W1y,H1y]=normalize(W1y,H1y);
        [W2y,H2y]=normalize(W2y,H2y);
        [W1,H1]=normalize(W1,H1);
        [W2,H2]=normalize(W2,H2);
        Wy=face_split(W1y,W2y);
        Hy=face_split(H1y,H2y);
        A=Hy'*Hy;
        B=X*Hy;
        lW=norm(A,2);
        alpha=tau/lW;
        
        % W-update
        for k=1:Iter_W
            G=Wy*A-B; 
            [W1y,W2y]=Upd_manifold(W1y,W2y,G,alpha);
            Wy=face_split(W1y,W2y);
        end
        W1n=W1y; W2n=W2y;
        W1y=W1n+beta*(W1n-W1);
        W2y=W2n+beta*(W2n-W2);

        [H1y,W1y]=normalize(H1y,W1y);
        [H2y,W2y]=normalize(H2y,W2y);
        W1n=W1n*diag(1./vecnorm(W1n));
        W2n=W2n*diag(1./vecnorm(W2n));
        [H1,W1]=normalize(H1,W1);
        [H2,W2]=normalize(H2,W2);
        Wy=face_split(W1y,W2y);
        Hy=face_split(H1y,H2y);
        C=Wy'*Wy;
        D=X'*Wy;
        lH=norm(C,2);
        alpha=tau/lH;

        % H-update
        for k=1:Iter_H
            G=Hy*C-D;
            [H1y,H2y]=Upd_manifold(H1y,H2y,G,alpha);
            Hy=face_split(H1y,H2y);
        end
        H1n=H1y; H2n=H2y;
        H1y=H1n+beta*(H1n-H1);
        H2y=H2n+beta*(H2n-H2);
        Hy=face_split(H1y,H2y);

        % Store time and error
        j=j+1;
        time(j)=toc(t_iter)+time(j-1);
        Wn=face_split(W1n,W2n);
        err(j)=norm(X-Wn*Hy','fro');

        % Extrapolation
        if err(j)<err(j-1)
            betaold=beta;
            beta=min(betat,gamma*beta);
            betat=min(1,gammat*betat);
            W1=W1n; W2=W2n; H1=H1n; H2=H2n;
        else
            err(j)=err(j-1);
            betat=betaold;
            beta=beta/eta;
            W1y=W1; W2y=W2; H1y=H1; H2y=H2;
        end 
    end

    % Output
    info=struct('err',err(1:j),'time',time(1:j));

end