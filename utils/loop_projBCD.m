function [W1,H1,W2,H2,info]=loop_projBCD(X,W1,H1,W2,H2,opts)
%% loop_projBCD: Main cycle for HadDec - projBCD case
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
    gammat=extrapar(4); eta=extrapar(5); kappa=extrapar(6); betaold=beta;
    tau=opts.tau;
    Iter_W=opts.Iter_W;
    Iter_H=opts.Iter_H;
    j=1;
    time=zeros(maxit,1);
    init_time=toc(t_start);
    time(j)=init_time;
    maxtime=opts.maxtime-init_time;
    sparseflag=opts.sparsity;
    extrapolflag=0;

    % Relative error function for unit Frobenius norm X
    err=zeros(maxit,1);
    W=face_split(W1,W2);
    H=face_split(H1,H2);
    if sparseflag
        A=H'*H;
        B=X*H;
        C=W'*W;
        err(j)=sqrt(1-2*W(:)'*B(:)+C(:)'*A(:));
    else
        err(j)=norm(X-W*H','fro');
    end

    % Update function
    if opts.noloops
        proj_Bmr=@(A) projBmr_noloop(A);
    else
        proj_Bmr=@(A) projBmr_loop(A);
    end

    W1y=W1;
    W2y=W2;
    H1y=H1;
    H2y=H2;

    % Main cycle
    if sparseflag
        % Loop with computation of the error for a sparse matrix
        while j<maxit && err(j)>tol && time(j)<maxtime
            
            t_iter=tic;
            [W1y,H1p,vecnormH1y]=normalize(W1y,H1y);
            [W2y,H2p,vecnormH2y]=normalize(W2y,H2y);
            Hp=face_split(H1p,H2p);
            A=Hp'*Hp;
            B=X*Hp;
            lW=norm(A,2);
            stepw=tau/lW;
            
            % W-update
            for k=1:Iter_W
                Wy=face_split(W1y,W2y);
                G=Wy*A-B; 
                [W1y,W2y]=proj_Bmr(Wy-stepw*G);
            end
            W1n=W1y*diag(1./vecnormH1y);
            W2n=W2y*diag(1./vecnormH2y);
            W1y=W1n+beta*(W1n-W1);
            W2y=W2n+beta*(W2n-W2);
    
            [H1y,W1p,vecnormW1y]=normalize(H1y,W1y);
            [H2y,W2p,vecnormW2y]=normalize(H2y,W2y);
            Wp=face_split(W1p,W2p);
            C=Wp'*Wp;
            D=X'*Wp;
            lH=norm(C,2);
            steph=tau/lH;
    
            % H-update
            for k=1:Iter_H
                Hy=face_split(H1y,H2y);
                G=Hy*C-D;
                [H1y,H2y]=proj_Bmr(Hy-steph*G);
            end
            H1n=H1y*diag(1./vecnormW1y);
            H2n=H2y*diag(1./vecnormW2y);
            H1y=H1n+beta*(H1n-H1);
            H2y=H2n+beta*(H2n-H2);
            Wy=Wp.*reshape(vecnormW2y*vecnormW1y',1,[]);
            C=Wy'*Wy;
            D=X'*Wy;
    
            % Store time and error
            j=j+1;
            time(j)=toc(t_iter)+time(j-1);
            Hn=face_split(H1n,H2n);
            A=Hn'*Hn;
            err(j)=sqrt(1-2*Hn(:)'*D(:)+C(:)'*A(:));
    
            % Extrapolation
            if err(j)<err(j-1)
                betaold=beta;
                beta=min(betat,gamma*beta);
                betat=min(kappa^extrapolflag,gammat*betat);
                W1=W1n; W2=W2n; H1=H1n; H2=H2n;
            else
                err(j)=err(j-1);
                betat=betaold;
                beta=beta/eta;
                extrapolflag=extrapolflag+1;
                if extrapolflag>10000
                    beta=0;
                end
                W1y=W1; W2y=W2; H1y=H1; H2y=H2;
            end 
        end
    else 
        % Loop with a more accurate computation of the error 
        while j<maxit && err(j)>tol && time(j)<maxtime
    
            t_iter=tic;
            [W1y,H1p,vecnormH1y]=normalize(W1y,H1y);
            [W2y,H2p,vecnormH2y]=normalize(W2y,H2y);
            Hp=face_split(H1p,H2p);
            A=Hp'*Hp;
            B=X*Hp;
            lW=norm(A,2);
            stepw=tau/lW;
            
            % W-update
            for k=1:Iter_W
                Wy=face_split(W1y,W2y);
                G=Wy*A-B; 
                [W1y,W2y]=proj_Bmr(Wy-stepw*G);
            end
            W1n=W1y*diag(1./vecnormH1y);
            W2n=W2y*diag(1./vecnormH2y);
            W1y=W1n+beta*(W1n-W1);
            W2y=W2n+beta*(W2n-W2);
    
            [H1y,W1p,vecnormW1y]=normalize(H1y,W1y);
            [H2y,W2p,vecnormW2y]=normalize(H2y,W2y);
            Wp=face_split(W1p,W2p);
            C=Wp'*Wp;
            D=X'*Wp;
            lH=norm(C,2);
            steph=tau/lH;
    
            % H-update
            for k=1:Iter_H
                Hy=face_split(H1y,H2y);
                G=Hy*C-D;
                [H1y,H2y]=proj_Bmr(Hy-steph*G);
            end
            H1n=H1y*diag(1./vecnormW1y);
            H2n=H2y*diag(1./vecnormW2y);
            H1y=H1n+beta*(H1n-H1);
            H2y=H2n+beta*(H2n-H2);
            Wy=Wp.*reshape(vecnormW2y*vecnormW1y',1,[]);
    
            % Store time and error
            j=j+1;
            time(j)=toc(t_iter)+time(j-1);
            Hn=face_split(H1n,H2n);
            err(j)=norm(X-Wy*Hn','fro');
    
            % Extrapolation
            if err(j)<err(j-1)
                betaold=beta;
                beta=min(betat,gamma*beta);
                betat=min(kappa^extrapolflag,gammat*betat);
                W1=W1n; W2=W2n; H1=H1n; H2=H2n;
            else
                err(j)=err(j-1);
                betat=betaold;
                beta=beta/eta;
                if extrapolflag>10000
                    beta=0;
                end
                W1y=W1; W2y=W2; H1y=H1; H2y=H2;
            end 
        end
    end  

    % Output
    info=struct('err',err(1:j),'time',time(1:j));

end