 
    % Dimensions of the problem
    n=4;
    m=n;
    r=floor(sqrt(min(n,m)));
    eta=1e-5;
    nmax=4;

%%
    % Example type
    example='oneseye';
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

    opts=struct('rank',r,'maxit',1,'init','FS3','tol',1e-16);
    warning('off', 'manopt:getHessian:approx') 
    tic;
    [X1,X2,err]=Had_manopt2(X,opts);
    timeHad=toc;
    finerr=0.5*norm(X-prodsvd(X1).*prodsvd(X2),'fro')^2;

    lw=1.3;
    semilogy(err,'r-o','LineWidth',lw)

function X=prodsvd(A)
    
    X=A.U*A.S*A.V';

end