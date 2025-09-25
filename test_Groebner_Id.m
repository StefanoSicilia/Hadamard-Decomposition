    
    r=3;
    n=r^2;
    I=eye(r);
    J=flip(I);
    one=ones(r,1);
    F=ones(r);
    C=diag(ones(r-1,1),1); C(r,1)=1;
    % r even case
    % X=kron(one,J);
    % U=[J;eye(r)];

    W=kron(one*one',I);
    H1=kron(F-I,J)+kron(I,I);
    H2=zeros(r^2);
    for i=1:r
        for j=1:i
            H2((i-1)*r+1:i*r,(j-1)*r+1:j*r)=C^(i-j);
        end
    end
    H2=H2+H2'-kron(I,I);
    H3=zeros(r^2);
    for i=1:r
        for j=1:i
            H3((i-1)*r+1:i*r,(j-1)*r+1:j*r)=C^(j-i);
        end
    end
    H3=H3+H3'-kron(I,I);
    A=W.*H2;




    % X=sym("X",[n,r]);
    % U=sym("U",[n,r]);
    % B=(X*X.').*(U*U.');
    % A=eye(n);
    % B=triu(B-A);
    % B=B(:);
    % a=gbasis(B.',[X(:),U(:)]);