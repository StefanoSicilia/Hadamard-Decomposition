function [W1,H1,W2,H2,info]=loop_proj(X,W1,H1,W2,H2,opts)
%% loop_proj: Main cycle for HadDec - proj case
% Performs Riemannian gradient block coordinate descent for Hadamard 
% decomposition of X, with rank provided by the starting approximation 
% X~=(W1*H1').*(W2*H2'). It implements alternative projections onto Bmr.
% See HadDec for more details.

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

    % Main cycle
    while j<maxit && err(j)>tol && time(j)<maxtime

        t_iter=tic;
        [W1,H1]=normalize(W1,H1);
        [W2,H2]=normalize(W2,H2);
        W=face_split(W1,W2);
        H=face_split(H1,H2);
        A=H'*H;
        B=X*H;
        lW=norm(A,2);
        alpha=tau/lW;
        
        % W-update
        for k=1:Iter_W
            W1old=W1;
            W2old=W2;
            G=W*A-B;
            [W1,W2]=proj_Bmr(W-alpha*G);
            W=face_split(W1,W2);
        end
        W1=W1+beta*(W1-W1old);
        W2=W2+beta*(W2-W2old);

        [H1,W1]=normalize(H1,W1);
        [H2,W2]=normalize(H2,W2);
        W=face_split(W1,W2);
        H=face_split(H1,H2);
        C=W'*W;
        D=X'*W;
        lH=norm(C,2);
        alpha=tau/lH;

        % H-update
        for k=1:Iter_H
            H1old=H1;
            H2old=H2;
            G=H*C-D;
            [H1,H2]=proj_Bmr(H-alpha*G);
            H=face_split(H1,H2); 
        end
        H1=H1+beta*(H1-H1old);
        H2=H2+beta*(H2-H2old);
        H=face_split(H1,H2);

        % Store time and error
        j=j+1;
        time(j)=toc(t_iter)+time(j-1);
        err(j)=norm(X-W*H','fro');

        % Extrapolation
        if err(j)<err(j-1)
            betaold=beta;
            beta=min(betat,gamma*beta);
            betat=min(1,gammat*betat);
        else
            betat=betaold;
            beta=beta/eta;
        end 
    end

    % Output 
    info=struct('err',err(1:j),'time',time(1:j));

end