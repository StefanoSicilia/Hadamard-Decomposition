function [W1,H1,W2,H2,info]=loop_BCD(X,W1,H1,W2,H2,opts)
%% loop_BCD: Main cycle for HadDec - BCD case
% Performs Riemannian gradient descent for Hadamard decomposition. See
% HadDec for more details.

    % Initializations and parameters
    tstart=tic;
    maxit=opts.maxit;
    extrapar=opts.momentum;
    beta=extrapar(1); betat=extrapar(2); gamma=extrapar(3); 
    gammat=extrapar(4); eta=extrapar(5); betaold=beta;
    err=zeros(maxit,1);
    time=zeros(maxit,1);
    init_time=toc(tstart);
    err(1)=norm(X-(W1*H1').*(W2*H2'),'fro');
    j=1;
    time(1)=init_time;
    maxtime=opts.maxtime-init_time;
    toler=opts.toler;
    reldiff(1:10)=toler+1;

    % If not, make X a full matrix
    if issparse(X)
        X=full(X);
    end

    % Update function
    if opts.noloops
        Upd_Fact=@(X,W1,H1,W2) UpdFact_noloop(X,W1,H1,W2);
    else
        Upd_Fact=@(X,W1,H1,W2) UpdFact_loop(X,W1,H1,W2);
    end

    % Main cycle
    while j<=maxit && err(j)>opts.tol && reldiff(j)>toler && time(j)<maxtime
        j=j+1;
        t_iter=tic;
        H2old=H2; H1old=H1; W2old=W2; W1old=W1;
        H2=Upd_Fact(X,W1,H1,W2);
        H2=H2+beta*(H2-H2old);
        W2=Upd_Fact(X',H1,W1,H2);
        W2=W2+beta*(W2-W2old);
        H1=Upd_Fact(X,W2,H2,W1);
        H1=H1+beta*(H1-H1old);
        W1=Upd_Fact(X',H2,W2,H1);
        W1=W1+beta*(W1-W1old);
        time(j)=toc(t_iter)+time(j-1);
        err(j)=norm(X-(W1*H1').*(W2*H2'),'fro');
        if err(j)<err(j-1)
            betaold=beta;
            beta=min(betat,gamma*beta);
            betat=min(1,gammat*betat);
        else
            betat=betaold;
            beta=beta/eta;
            H2=H2old; H1=H1old; W2=W2old; W1=W1old;
        end
        if j>10
            reldiff(j)=abs(err(j-10)-err(j))/err(j-10);
        end
    end

    % Output
    info=struct('err',err(1:j),'time',time(1:j));

end