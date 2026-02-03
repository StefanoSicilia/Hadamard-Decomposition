function [W1,H1,W2,H2,info]=loop_BCD(X,W1,H1,W2,H2,opts)
%% loop_BCD: Main cycle for HadDec
% Performs Riemannian gradient descent for Hadamard decomposition. See
% HadDec_BCD for more details.

    % Initializations and parameters
    tstart=tic;
    maxit=opts.maxit;
    extrapar=opts.momentum;
    beta=extrapar(1); betat=extrapar(2); gamma=extrapar(3); 
    gammat=extrapar(4); eta=extrapar(5);
    err=zeros(maxit,1);
    time=zeros(maxit,1);
    err(1)=norm(X-(W1*H1').*(W2*H2'),'fro');
    i=1;
    init_time=toc(tstart);
    time(1)=init_time;
    maxtime=opts.maxtime-init_time;
    tstart=tic;

    % Main cycle
    while i<=maxit && err(i)>opts.tol && time(i)<maxtime
        i=i+1;
        H2old=H2; H1old=H1; W2old=W2; W1old=W1;
        H2=UpdFact(X,W1,H1,W2,H2);
        H2=H2+beta*(H2-H2old);
        W2=UpdFact(X',H1,W1,H2,W2);
        W2=W2+beta*(W2-W2old);
        H1=UpdFact(X,W2,H2,W1,H1);
        H1=H1+beta*(H1-H1old);
        W1=UpdFact(X',H2,W2,H1,W1);
        W1=W1+beta*(W1-W1old);
        err(i)=norm(X-(W1*H1').*(W2*H2'),'fro');
        time(i)=toc(tstart);
        if err(i)<err(i-1)
            beta=min(betat,gamma*beta);
            betat=min(1,gammat*betat);
        else
            beta=beta/eta;
            betat=beta;
        end
    end

    % Output
    err(i)=norm(X-(W1*H1').*(W2*H2'),'fro');
    info=struct('err',err(1:i),'time',time(1:i));


end
