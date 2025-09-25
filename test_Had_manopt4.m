%% Not working :(

    % Dimensions of the problem
    n=4;
    m=5;
    r=floor(sqrt(min(n,m)));
    eta=1e-5;
    nmax=4;

%%
    % Example type
    example='eye';
    rng(1)
    imageflag=0;
    switch example
        case 'eye'
            X=eye(n,m);
        case 'oneseye'
            X=ones(n,m)-eye(n,m);
        case 'random'
            X=rand(n,m);
        case 'randexact'
            W1=randi(nmax,n,r);
            H1=randi(nmax,m,r);
            W2=randi(nmax,n,r);
            H2=randi(nmax,m,r);
            X=(W1*H1').*(W2*H2');
        case 'randexactpert'
            W1=randi(nmax,n,r);
            H1=randi(nmax,m,r);
            W2=randi(nmax,n,r);
            H2=randi(nmax,m,r);
            X=(W1*H1').*(W2*H2')+eta*randi(nmax,n,m);
        case 'football'
            n=115;
            m=n;
            r=9;
            X=full(load("./datasets/football.mat").A);
        case 'cat'
            n=400;
            m=600;
            r=20;
            X=double(imread('./datasets/cat.jpg'));
            X=X+1*randn(n,m);
            imageflag=1;
        case 'cameraman'
            n=225;
            m=n;
            r=15;
            X=double(imread('./datasets/cameraman.jpg'));
            X=X(:,:,1);
            %A=A+1*randn(n,m);
            imageflag=1;
        case 'peppers'
            n=192;
            m=204;
            r=13;
            X=double(imread('./datasets/peppers.jfif'));
            X=255-X(:,:,1);
            imageflag=1;
        case 'apple'
            n=180;
            m=281;
            r=13;
            X=double(imread('./datasets/apple.jfif'));
            X=X(:,:,1);
            imageflag=1;
        case 'blender'
            n=512;
            m=n;
            r=22;
            X=double(imread('./datasets/blender.png'));
            imageflag=1;
        otherwise
            error('Example type not available.')
    end
    normX=norm(X,'fro');
    %imshow(255*X/normX);
    X=X/normX;

    opts=struct('rank',r,'init','FS2','maxtime',10);
    warning('off', 'manopt:getHessian:approx') 
    tic;
    [X1,X2,err]=Had_manopt4(X,opts);
    timeHad=toc;
    finerr=0.5*norm(X-prodsvd(X1).*prodsvd(X2),'fro')^2;

    if imageflag
        sq=255;
        figure(1)
        imshow(sq*prodsvd(X1).*prodsvd(X2))
        figure(2)
        [U,S,V]=svd(A);
        r2=2*r;
        Sr=sqrt(S(1:r2,1:r2));
        U=U(:,1:r2)*Sr;
        V=V(:,1:r2)*Sr;
        imshow(255*U*V')
        fin_errSVDmatlab=0.5*norm(A-U*V','fro')^2;
        figure(3)
        imshow(255*A)
        s1=ssim(A,(X_svd*Y_svd').*(U_svd*V_svd'));
        sW1=ssim(A,(X_Wsvd*Y_Wsvd').*(U_Wsvd*V_Wsvd'));
    end

function X=prodsvd(A)
    
    X=A.U*A.S*A.V';

end