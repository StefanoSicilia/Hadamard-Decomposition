%% Script to test things on Hadamard product

    r=2;
    n=r^2;
    m=n;
    nmax=2;
    rng(1)
    epsilon=1;
    XE=randi(nmax,n,r);
    XE(2:end,:)=zeros(n-1,r);
    rng(21)
    X=randi(nmax,n,r)+epsilon*XE;
    U=randi(nmax,n,r);
    Y=randi(nmax,m,r);
    V=randi(nmax,m,r);
    A=(X*Y').*(U*V');

    Z=column_Had(X,U);
    Q=column_Had(Y,V);
    err=norm(A-Z*Q','fro');
    checkrow=det(reshape(Z(1,:),r,r));
    normcheckrow=norm(checkrow);
    checkcol=det(reshape(Q(1,:),r,r));
    normcheckcol=norm(checkcol);

    S=randi(nmax,r,r);
    Sigma=randi(nmax,r,r);
    A1=(X*S*Y').*(U*Sigma*V');
    A2=face_split(X,U)*kron(S,Sigma)*face_split(Y,V)';

    % [X,~]=qr(X,0);
    % [U,~]=qr(U,0);
    % G=face_split(X,U);

    % B=zeros(n,m);
    % for i=1:n
    %     for j=1:m
    %         B(i,j)=norm(X(i,:))*norm(Y(j,:))*norm(U(i,:))*norm(V(j,:));
    %     end
    % end
    % erra='------';
    % errb='------';
    % i=1; j=3; h=4; k=2;
    % xi=X(i,:); yj=Y(j,:); uh=U(h,:); vk=V(k,:);
    % xi=xi/norm(xi);
    % yj=yj/norm(yj);
    % uh=uh/norm(uh);
    % vk=vk/norm(vk);
    % %gyj=yj;
    % gyj=yj-yj*xi'/norm(xi)^2*xi; gyj=gyj/norm(gyj);
    % %puh=uh;
    % %pvk=vk;
    % puh=(uh*xi')*xi+(uh*gyj')*gyj; puh=puh/norm(puh);
    % pvk=(vk*xi')*xi+(vk*gyj')*gyj; pvk=pvk/norm(pvk);
    % thetaij=xi*yj';
    % thetahk=puh*pvk';
    % thetaik=xi*pvk';
    % thetahj=puh*yj';
    % errangle1=1-cos(acos(thetaik)+acos(thetahj)-acos(thetaij)-acos(thetahk));
    % errangle2=1-cos(acos(thetaik)-acos(thetahj)+acos(thetaij)-acos(thetahk));
    % errangle3=1-cos(acos(thetaik)-acos(thetahj)-acos(thetaij)+acos(thetahk));
    % errangle4=1-cos(acos(thetaik)+acos(thetahj)-acos(thetaij)+acos(thetahk));
    % errangle5=1-cos(-acos(thetaik)+acos(thetahj)+acos(thetaij)+acos(thetahk));
    % errangle6=1-cos(acos(thetaik)-acos(thetahj)+acos(thetaij)+acos(thetahk));
    % errangle7=1-cos(acos(thetaik)+acos(thetahj)+acos(thetaij)-acos(thetahk));
    % errangle=1-cos(acos(B(i,k))+acos(B(h,j))-acos(B(i,j))-acos(B(h,k)));

    rng(2)
    E=randi(nmax,n,r);
    G=randi(nmax,n,r);
    F=randi(nmax,m,r);
    H=randi(nmax,m,r);
    C=face_split(Y,V)'*face_split(F,H)
    D=kron(Y'*H,V'*F);


       
    
