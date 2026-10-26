function [W1,H1,W2,H2,info]=loop_Manopt_WH(X,W1,H1,W2,H2,opts)
%% loop_Manopt_WH: Main cycle for HadDec - Manopt_WH case
% Performs Riemannian gradient descent for Hadamard decomposition on W and
% H constrained to be on the manifolds Bmr and Bnr. See HadDec for more
% details.

    % Initialization and parameters
    tstart=tic;
    [m,n]=size(X);
    r1=size(W1,2);
    r2=size(W2,2);

    M1=fixedrankembeddedfactory(r2,r1,1);
    MW=powermanifold(M1,m);
    MH=powermanifold(M1,n);
    M=productmanifold(struct('W',MW,'H',MH));
    problem.M=M;

    P0.W=make_starting_point(W1,W2,M1);
    P0.H=make_starting_point(H1,H2,M1);
    toler=opts.toler;
    options.verbosity=0;
    options.tolcost=0.5*(opts.tol)^2;
    options.maxiter=Inf;
    options.tolgradnorm=0.5*(opts.tol)^2;
    init_time=toc(tstart);
    options.maxtime=opts.maxtime-init_time;
    options.stopfun=@(problem,x,info,last) reldiffStop(problem,x,info,last,toler);
    function stopnow=reldiffStop(~,~,info,last,toler)
        delay=10;
        startcheck=50;
        if last>=startcheck
            err_curr=sqrt(2*info(last).cost);
            err_prev=sqrt(2*info(last-delay).cost);
            stopnow=abs(err_prev-err_curr)/err_prev<toler;
        else
            stopnow=false;
        end
    end

    % Calling Manopt's trustregions. For sparse X, use matrix-free
    % expressions so that the dense product W*H' and residual W*H'-X are
    % never formed.
    if issparse(X)
        Xt=X';
        Xnorm2=norm(X,'fro')^2;
        problem.cost=@(P)costfun_sparse(P,X,M1,Xnorm2);
        problem.egrad=@(P)egradfun_sparse(P,X,Xt,M1);
        problem.ehess=@(P,eta)ehessfun_sparse(P,eta,X,Xt,M1);
    else
        problem.cost=@(P)costfun(P,X,M1);
        problem.egrad=@(P)egradfun(P,X,M1);
        problem.ehess=@(P,eta)ehessfun(P,eta,X,M1);
    end
    % warning('off', 'manopt:getHessian:approx')
    warning('on')
    [P,~,info] = trustregions(problem,P0,options);

    % Output
    info=struct('err',sqrt(2*[info.cost]'),'time',[info.time]');
    [W1,W2]=extract_factors(P.W);
    [H1,H2]=extract_factors(P.H);

end



function C=make_starting_point(A1,A2,M1)
    p=size(A1,1);
    C=cell(p,1);
    for i=1:p
        A=A2(i,:)'*A1(i,:);
        C{i}=M1.rand();
        [U,S,V]=svd(A,'econ');
        C{i}.U=U(:,1);
        C{i}.S=S(1,1);
        C{i}.V=V(:,1);
    end
end

function f=costfun(P,X,M1)
    W=points_to_matrix(P.W,M1);
    H=points_to_matrix(P.H,M1);
    E=W*H'-X;
    f=0.5*norm(E,'fro')^2;
end

function g=egradfun(P,X,M1)
    W=points_to_matrix(P.W,M1);
    H=points_to_matrix(P.H,M1);
    E=W*H'-X;
    GW=E*H;
    GH=E'*W;
    gW=cell(size(P.W));
    gH=cell(size(P.H));
    for i=1:length(P.W)
        gW{i}=reshape(GW(i,:),size(P.W{i}.U,1),size(P.W{i}.V,1));
    end
    for i=1:length(P.H)
        gH{i}=reshape(GH(i,:),size(P.H{i}.U,1),size(P.H{i}.V,1));
    end
    g=struct('W',{gW},'H',{gH});
end

function h=ehessfun(P,eta,X,M1)
    W=points_to_matrix(P.W,M1);
    H=points_to_matrix(P.H,M1);
    dW=tangents_to_matrix(P.W,eta.W,M1);
    dH=tangents_to_matrix(P.H,eta.H,M1);
    E=W*H'-X;
    % dGW=dW*(H'*H)+E*dH;
    % dGH=dH*(W'*W)+E'*dW;

    dE = dW*H' + W*dH';
    dGW = dE*H + E*dH;
    dGH = dE'*W + E'*dW;
    hW=cell(size(P.W));
    hH=cell(size(P.H));
    for i=1:length(P.W)
        hW{i}=reshape(dGW(i,:),size(P.W{i}.U,1),size(P.W{i}.V,1));
    end
    for i=1:length(P.H)
        hH{i}=reshape(dGH(i,:),size(P.H{i}.U,1),size(P.H{i}.V,1));
    end
    h=struct('W',{hW},'H',{hH});
end

function f=costfun_sparse(P,X,M1,Xnorm2)
    W=points_to_matrix(P.W,M1);
    H=points_to_matrix(P.H,M1);

    % ||X-WH'||_F^2 = ||X||_F^2 - 2<W,XH>
    %                    + <W'W,H'H>.
    % X*H uses sparse matrix multiplication and WH' is never formed.
    XH=X*H;
    WtW=W'*W;
    HtH=H'*H;
    residual_norm2=Xnorm2-2*sum(W.*XH,'all')+sum(WtW.*HtH,'all');

    % Guard against a tiny negative value caused by roundoff cancellation.
    f=0.5*max(real(residual_norm2),0);
end

function g=egradfun_sparse(P,X,Xt,M1)
    W=points_to_matrix(P.W,M1);
    H=points_to_matrix(P.H,M1);

    % Joint Euclidean gradient, evaluated without forming WH'.
    GW=W*(H'*H)-X*H;
    GH=H*(W'*W)-Xt*W;
    gW=cell(size(P.W));
    gH=cell(size(P.H));
    for i=1:length(P.W)
        gW{i}=reshape(GW(i,:),size(P.W{i}.U,1),size(P.W{i}.V,1));
    end
    for i=1:length(P.H)
        gH{i}=reshape(GH(i,:),size(P.H{i}.U,1),size(P.H{i}.V,1));
    end
    g=struct('W',{gW},'H',{gH});
end

function h=ehessfun_sparse(P,eta,X,Xt,M1)
    W=points_to_matrix(P.W,M1);
    H=points_to_matrix(P.H,M1);
    dW=tangents_to_matrix(P.W,eta.W,M1);
    dH=tangents_to_matrix(P.H,eta.H,M1);

    % Directional derivative of the joint Euclidean gradient. These
    % expressions use sparse products with X and never form WH' or d(WH').
    HtH=H'*H;
    WtW=W'*W;
    dHtH=dH'*H+H'*dH;
    dWtW=dW'*W+W'*dW;
    dGW=dW*HtH+W*dHtH-X*dH;
    dGH=dH*WtW+H*dWtW-Xt*dW;
    hW=cell(size(P.W));
    hH=cell(size(P.H));
    for i=1:length(P.W)
        hW{i}=reshape(dGW(i,:),size(P.W{i}.U,1),size(P.W{i}.V,1));
    end
    for i=1:length(P.H)
        hH{i}=reshape(dGH(i,:),size(P.H{i}.U,1),size(P.H{i}.V,1));
    end
    h=struct('W',{hW},'H',{hH});
end

function A=points_to_matrix(C,M1)
    p=length(C);
    A=zeros(p,size(C{1}.U,1)*size(C{1}.V,1));
    for i=1:p
        Ai=M1.triplet2matrix(C{i});
        A(i,:)=Ai(:)';
    end
end

function A=tangents_to_matrix(X,Z,M1)
    p=length(Z);
    A=zeros(p,size(X{1}.U,1)*size(X{1}.V,1));
    for i=1:p
        Zi=M1.tangent2ambient(X{i},Z{i});
        Zi=M1.triplet2matrix(Zi);
        A(i,:)=Zi(:)';
    end
end

function [A1,A2]=extract_factors(C)
    p=length(C);
    r1=size(C{1}.V,1);
    r2=size(C{1}.U,1);
    A1=zeros(p,r1);
    A2=zeros(p,r2);
    for i=1:p
        U=C{i}.U;
        S=C{i}.S;
        V=C{i}.V;
        A1(i,:)=(S(1,1)*V(:,1))';
        A2(i,:)=U(:,1)';
    end
end
