function [W1,H1,W2,H2,info]=loop_manBCD_sparse(X,W1,H1,W2,H2,opts)
%% loop_manBCD_sparse: Main cycle for HadDec - manBCDsparse case
% Performs Riemannian gradient descent for Hadamard decomposition using 
% sparse operations. See HadDec for more details.

    % Initialization and parameters
    t_start=tic;
    maxit=opts.maxit;
    tol=opts.tol;
    extrapar=opts.momentum;
    beta=extrapar(1); betat=extrapar(2); gamma=extrapar(3); 
    gammat=extrapar(4); eta=extrapar(5);
    tau=opts.tau;

    W=face_split(W1,W2);
    H=face_split(H1,H2);

    Hblock=opts.Hblock;
    Wblock=opts.Wblock;
    err=zeros(maxit,1);
    XH=X*H;
    WtW=W'*W;
    HtH=H'*H;
    err(1)=0.5*(1-2*trace(XH'*W)+trace(WtW*HtH));
    i=1;
    j=1;
    time=zeros(maxit,1);
    init_time=toc(t_start);
    time(1)=init_time;
    maxtime=opts.maxtime-init_time;
    

    % Main cycle
    while j<maxit && err(i)>tol && time(i)<maxtime
        % W-update
        for k=1:Wblock
            t_iter=tic;
            W1old=W1;
            W2old=W2;
            G=W*HtH-X*H;
            alpha=tau/norm(H'*H,2); 
            [W1,W2]=Upd_manifold(W1,W2,G,alpha);
            W1=W1+beta*(W1-W1old);
            W2=W2+beta*(W2-W2old);
            W=face_split(W1,W2);
            i=i+1;
            XtW=X'*W;
            WtW=W'*W;
            time(i)=toc(t_iter)+time(i-1);
            err(i)=0.5*(1-2*trace(H'*XtW)+trace(WtW*HtH));
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
            t_iter=tic;
            H1old=H1;
            H2old=H2;
            G=H*WtW-XtW;
            alpha=tau/norm(W'*W,2);
            [H1,H2]=Upd_manifold(H1,H2,G,alpha);
            H1=H1+beta*(H1-H1old);
            H2=H2+beta*(H2-H2old);
            H=face_split(H1,H2);
            i=i+1;
            XH=X*H;
            HtH=H'*H;
            time(i)=toc(t_iter)+time(i-1);
            err(i)=0.5*(1-2*trace(XH'*W)+trace(WtW*HtH));
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
    err(i)=0.5*(1-2*trace(XH'*W)+trace(WtW*HtH));
    info=struct('err',sqrt(2*err(1:i)),'time',time(1:i));

end