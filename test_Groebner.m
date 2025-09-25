
    % n=4;
    % m=4;
    % r=2;
    % X=sym("X",[n,r]);
    % Y=sym("Y",[m,r]);
    % U=sym("U",[n,r]);
    % V=sym("V",[m,r]);
    % B=tril((X*X.').*(U*U.'),0);
    % %A=ones(n,m)-eye(n,m);
    % %A=ones(n,m);
    % %A=sym("A",[n,m]); A = tril(A,0);
    % A=ones(n,m); A = tril(A,0);
    % B=B(:);
    % %a=gbasis(B.',[X(:);U(:)]);
    % solution=solve(B==A(:),[X(:);U(:)]);

    % n=2;
    % m=2;
    % r=1;
    % X=sym("X",[n,r]);
    % Y=sym("Y",[m,r]);
    % U=sym("U",[n,r]);
    % V=sym("V",[m,r]);
    % B=(X*Y.').*(U*V.');
    % A=ones(n,m);
    % %A=sym("A",[n,m]);
    % B=B-A;
    % B=B(:);
    % eqn=B==0;
    % a=solve(eqn,[X(:);Y(:);U(:);V(:)]);

    % n=4;
    % m=4;
    % r=2;
    % xi=sym("xi",[r,1]);
    % ui=sym("ui",[r,1]);
    % ai=sym("ai",[r,1]);
    % bi=sym("bi",[r,1]);
    % C=sym("C",[m,r]);
    % D=sym("D",[m,r]);
    % system=sym("system",[n,1]);
    % alpha=sym("alpha",[n,1]);
    % for j=1:n
    %     cj=C(j,:).';
    %     dj=D(j,:).';
    %     system(j)=(xi.'*cj)*(ui.'*dj);
    %     alpha(j)=(ai.'*cj)+(bi.'*dj);
    % end
    % solution=solve(system==alpha,[xi,ui]);
    %grobner_basis=gbasis((system-alpha).',[xi,ui]);

    % n=4;
    % m=4;
    % r=2;
    % ai=sym("ai",[r,1]);
    % bi=sym("bi",[r,1]);
    % cj=sym("cj",[r,1]);
    % dj=sym("dj",[r,1]);
    % xi=ai./bi+bi./ai;
    % yj=cj.*cj+dj.*dj;
    % ui=ai.*ai-ai.*bi;
    % vj=cj.*cj-dj.*dj;
    % res=(xi.'*yj)*(ui.'*vj)-(ai.'*cj)-(bi.'*dj);

    % n=4;
    % m=4;
    % r=2;
    % W=sym("W",[n,m]);
    % H=sym("H",[n,m]);
    % system=sym("system",[16+8,1]);
    % for i=1:n
    %     system(i)=det(reshape(W(i,:),r,r));
    % end
    % for i=1:m
    %     system(i+n)=det(reshape(H(:,i),r,r));
    % end
    % system(n+m+1:end)=reshape(W*H.',n*m,1);
    % A=sym("A",[n,m]); a=[zeros(n+m,1);reshape(A,n*m,1)];
    % %solution=solve(system==a,[W(:),H(:)]);
    % grobner_basis=gbasis((system-a).',[W(:),H(:)]);

    % n=4;
    % m=4;
    % r=2;
    % P=sym("P",[n,m]);
    % Q=sym("Q",[n,m]);
    % 
    % P(1,2)=0; P(1,4)=0;
    % P(2,1)=0; P(2,3)=0;
    % P(3,3)=0; P(3,4)=0;
    % P(4,1)=0; P(4,2)=0;
    % 
    % Q(1,3)=0; Q(1,4)=0;
    % Q(2,1)=0; Q(2,2)=0;
    % Q(3,2)=0; Q(3,4)=0;
    % Q(4,1)=0; Q(4,3)=0;
    % 
    % system=reshape(P*Q.',n*m,1);
    % A=sym("A",[n,m]); 
    % A(1,4)=0; A(2,3)=0; A(3,2)=0; A(4,1)=0; 
    % a=reshape(A,n*m,1);
    % %A=ones(n,m); a=reshape(A,n*m,1);
    % %A=ones(n,m)-eye(n,m);
    % vecP=P(:);
    % p=[vecP(1);vecP(3);vecP(6:7);vecP(9);vecP(12);vecP(14);vecP(16)];
    % vecQ=Q(:);
    % q=[vecQ(1);vecQ(3);vecQ(5);vecQ(8);vecQ(10:11);vecQ(14);vecQ(16)];
    % solution=solve(system==a,[p,q],Real=true);
    % grobner_basis=gbasis((system-a).',[p,q]);

    % n=4;
    % m=4;
    % r=2;
    % P=sym("P",[n,m]);
    % Q=sym("Q",[n,m]);
    % 
    % P(1,2)=0; P(1,4)=0;
    % P(2,1)=0; P(2,3)=0;
    % P(3,3)=0; P(3,4)=0;
    % P(4,1)=0; P(4,2)=0;
    % 
    % Q(1,1)=0; Q(1,2)=0;
    % Q(2,3)=0; Q(2,4)=0;
    % Q(3,1)=0; Q(3,3)=0;
    % Q(4,2)=0; Q(4,4)=0;
    % 
    % system=reshape(P*Q.',n*m,1);
    % A=sym("A",[n,m]); 
    % A=randi(5,n,m);
    % A(1,3)=0; A(2,4)=0; A(3,1)=0; A(4,2)=0; 
    % a=reshape(A,n*m,1);
    % %A=ones(n,m); a=reshape(A,n*m,1);
    % %A=ones(n,m)-eye(n,m);
    % vecP=P(:);
    % p=[vecP(1);vecP(3);vecP(6:7);vecP(9);vecP(12);vecP(14);vecP(16)];
    % vecQ=Q(:);
    % q=[vecQ(2);vecQ(4);vecQ(6:7);vecQ(9);vecQ(12);vecQ(13);vecQ(15)];
    % solution=solve(system==a,[p,q],Real=true);
    % %grobner_basis=gbasis((system-a).',[p,q]);
    % if size(solution.P1_1,1)>1
    %     Pstar=[solution.P1_1(1),0,solution.P1_3(1),0;
    %         0,solution.P2_2(1),0,solution.P2_4(1);
    %         solution.P3_1(1),solution.P3_2(1),0,0;
    %         0,0,solution.P4_3(1),solution.P4_4(1);];
    %     Qstar=[0,0,solution.Q1_3(1),solution.Q1_4(1);
    %         solution.Q2_1(1),solution.Q2_2(1),0,0;
    %         0,solution.Q3_2(1),0,solution.Q3_4(1);
    %         solution.Q4_1(1),0,solution.Q4_3(1),0;];
    % else
    %     Pstar=[solution.P1_1,0,solution.P1_3,0;
    %         0,solution.P2_2,0,solution.P2_4;
    %         solution.P3_1,solution.P3_2,0,0;
    %         0,0,solution.P4_3,solution.P4_4;];
    %     Qstar=[0,0,solution.Q1_3,solution.Q1_4;
    %         solution.Q2_1,solution.Q2_2,0,0;
    %         0,solution.Q3_2,0,solution.Q3_4;
    %         solution.Q4_1,0,solution.Q4_3,0;];
    % end
    % Astar=Pstar*Qstar.';
    % norm(Astar-A,'fro')

    % n=4;
    % m=4;
    % r=2;
    % P=sym("P",[n,m]);
    % Q=sym("Q",[n,m]);
    % 
    % P(1,3)=0; P(1,4)=0;
    % P(2,3)=0; P(2,4)=0;
    % P(3,1)=0; P(3,2)=0;
    % P(4,1)=0; P(4,2)=0;
    % 
    % Q(1,1)=0; Q(1,3)=0;
    % Q(2,2)=0; Q(2,4)=0;
    % Q(3,1)=0; Q(3,3)=0;
    % Q(4,2)=0; Q(4,4)=0;
    % 
    % system=reshape(P*Q.',n*m,1);
    % A=sym("A",[n,m]); 
    % a=reshape(A,n*m,1);
    % %A=ones(n,m); a=reshape(A,n*m,1);
    % %A=ones(n,m)-eye(n,m);
    % vecP=P(:);
    % p=vecP(vecP~=0);
    % vecQ=Q(:);
    % q=vecQ(vecQ~=0);
    % solution=solve(system==a,[p,q],Real=true);
    % grobner_basis=gbasis((system-a).',[p,q]);

    % n=4;
    % m=4;
    % r=2;
    % P=sym("P",[n,r]);
    % Q=sym("Q",[n,r]);
    % lambda=sym("lambda",[2,1]);
    % mu=sym("mu",[2,1]);
    % 
    % I=eye(r);
    % O=zeros(r);
    % L=[lambda(1)*I,O;O,lambda(2)*I];
    % Pbig=[P,L*P];
    % 
    % M=[mu(1)*I,O;O,mu(2)*I];
    % Qbig=[Q,M*Q];
    % system=reshape(Pbig*Qbig.',n*m,1);
    % A=sym("A",[n,m]); A(1,2)=0; 
    % a=reshape(A,n*m,1);
    % %A=ones(n,m); a=reshape(A,n*m,1);
    % %A=ones(n,m)-eye(n,m);
    % vecP=P(:);
    % p=vecP(vecP~=0);
    % q=Q(:);
    % solution=solve(system==a,[p;q;lambda;mu],Real=true);
    % %grobner_basis=gbasis((system-a).',[p;q;lambda;mu]);
    % % solution=solve(system==a,[p;q;lambda],Real=true);
    % % grobner_basis=gbasis((system-a).',[p;q;lambda]);

    %% The case with 4 zeros
    % n=4;
    % m=4;
    % r=2;
    % P=sym("P",[n,m]);
    % Q=sym("Q",[n,m]);
    % 
    % P(1,3)=0; P(1,4)=0;
    % P(2,1)=0; P(2,2)=0;
    % P(3,2)=0; P(3,4)=0;
    % P(4,1)=0; P(4,3)=0;
    % 
    % Q(1,1)=0; Q(1,2)=0;
    % Q(2,3)=0; Q(2,4)=0;
    % Q(3,1)=0; Q(3,3)=0;
    % Q(4,2)=0; Q(4,4)=0;
    % 
    % system=reshape(P*Q.',n*m,1);
    % A=sym("A",[n,m]); 
    % A(1,1)=0; A(2,2)=0; A(3,3)=0; A(4,4)=0; 
    % a=reshape(A,n*m,1);
    % %A=ones(n,m); a=reshape(A,n*m,1);
    % A=ones(n,m)-eye(n,m);
    % vecP=P(:);
    % p=vecP(vecP~=0);
    % vecQ=Q(:);
    % q=vecQ(vecQ~=0);
    % solution=solve(system==a,[p,q],Real=true);
    % %grobner_basis=gbasis((system-a).',[p,q]);

    %% The case with 2 zeros
    n=4;
    m=4;
    r=2;
    X=sym("X",[n,r]);
    U=sym("U",[n,r]);
    Y=sym("Y",[m,r]);
    V=sym("V",[m,r]);

    X(3:4,:)=[1,0;0,1];
    Y(3:4,:)=[1,0;0,1];
    U(1:2,:)=[1,1;1,1];
    V(1:2,:)=[1,0;0,1]; 

    P=sym_fs(X,U);
    Q=sym_fs(Y,V);
    spy(P*Q.')

    system=reshape(P*Q.',n*m,1);
    A=sym("A",[n,m]); 
    A(1,2)=0; A(4,3)=0;
    a=reshape(A,n*m,1); 
    vecX=X(:);
    x=vecX(vecX~=real(vecX)); 
    vecU=U(:); 
    u=vecU(vecU~=real(vecU)); 
    vecY=Y(:);
    y=vecY(vecY~=real(vecY)); 
    vecV=V(:);
    v=vecV(vecV~=real(vecV)); 
    solution=solve(system==a,[x;u;y;v],Real=true);

    % Xstar=zeros(n,r); Xstar(3:4,:)=eye(2);
    % Ustar=zeros(n,r); Ustar(1:2,:)=ones(2);
    % Ystar=zeros(m,r); Ystar(3:4,:)=eye(2);
    % Vstar=zeros(m,r); Vstar(1:2,:)=eye(2); Vstar(3,2)=1; Vstar(4,2)=1;
    % rng(133)
    % A=rand(n,m); %A(3,4)=0; A(4,3)=0;
    % Xstar(1,1)=A(1,3)*dete(A,[3,4],[1,3])*dete(A,4,2)/(dete(A,[3,4],[1,2])*dete(A,4,3));
    % Xstar(1,2)=A(1,4)*dete(A,[3,4],[1,4])*dete(A,3,2)/(dete(A,[3,4],[1,2])*dete(A,3,4));
    % Xstar(2,1)=A(2,3)*dete(A,[3,4],[1,3])*dete(A,4,2)/(dete(A,[3,4],[1,2])*dete(A,4,3));
    % Xstar(2,2)=A(2,4)*dete(A,[3,4],[1,4])*dete(A,3,2)/(dete(A,[3,4],[1,2])*dete(A,3,4));
    % Ystar(1,1)=dete(A,[3,4],[2,3])*dete(A,4,3)/(dete(A,[3,4],[1,3])*dete(A,4,2));
    % Ystar(1,2)=-dete(A,[3,4],[2,4])*dete(A,3,4)/(dete(A,[3,4],[1,4])*dete(A,3,2));
    % Ystar(2,1)=dete(A,4,3)/dete(A,4,2);
    % Ystar(2,2)=-dete(A,3,4)/dete(A,3,2);
    % Ustar(3,1)=A(3,1)*dete(A,[3,4],[1,3])*dete(A,4,2)/(dete(A,[3,4],[2,3])*dete(A,4,3));
    % Ustar(4,1)=-A(4,1)*dete(A,[3,4],[1,4])*dete(A,3,2)/(dete(A,[3,4],[2,4])*dete(A,3,4));
    % Ustar(3,2)=A(3,2)*dete(A,4,2)/dete(A,4,3);
    % Ustar(4,2)=-A(4,2)*dete(A,3,2)/dete(A,3,4);
    % Vstar(3,1)=-dete(A,[3,4],[2,3])*dete(A,4,1)/(dete(A,[3,4],[1,3])*dete(A,4,2));
    % Vstar(4,1)=-dete(A,[3,4],[2,4])*dete(A,3,1)/(dete(A,[3,4],[1,4])*dete(A,3,2));
    % vecXstar=Xstar(:);
    % xstar=vecXstar(vecXstar~=0);
    % %subs(solution.V4_1,a,A(:))-Vstar(4,1)
    % norm(A-(Xstar*Ystar').*(Ustar*Vstar'))

    Xstar=zeros(n,r); Xstar(3:4,:)=eye(2);
    Ustar=zeros(n,r); Ustar(1:2,:)=ones(2);
    Ystar=zeros(m,r); Ystar(3:4,:)=eye(2);
    Vstar=zeros(m,r); Vstar(1:2,:)=eye(2); Vstar(3,2)=1; Vstar(4,2)=1;
    rng(133)
    A=rand(n,m); %A(3,4)=0; A(4,3)=0;
    Xstar(1,1)=A(1,3)*dete(A,[3,4],[1,3])/(dete(A,[3,4],[1,2]));
    Xstar(1,2)=A(1,4)*dete(A,[3,4],[1,4])/(dete(A,[3,4],[1,2]));
    Xstar(2,1)=A(2,3)*dete(A,[3,4],[1,3])/(dete(A,[3,4],[1,2]));
    Xstar(2,2)=A(2,4)*dete(A,[3,4],[1,4])/(dete(A,[3,4],[1,2]));
    Ystar(1,1)=dete(A,[3,4],[2,3])/(dete(A,[3,4],[1,3]));
    Ystar(1,2)=-dete(A,[3,4],[2,4])/(dete(A,[3,4],[1,4]));
    Ystar(2,1)=dete(A,4,3)/dete(A,4,2);
    Ystar(2,2)=-dete(A,3,4)/dete(A,3,2);
    Ystar(3,1)=dete(A,4,2)/dete(A,4,3);
    Ystar(4,2)=dete(A,3,2)/dete(A,3,4);
    Ustar(3,1)=A(3,1)*dete(A,[3,4],[1,3])/(dete(A,[3,4],[2,3]));
    Ustar(4,1)=-A(4,1)*dete(A,[3,4],[1,4])/(dete(A,[3,4],[2,4]));
    Ustar(3,2)=A(3,2);
    Ustar(4,2)=-A(4,2);
    Ustar(1,:)=dete(A,4,2)/dete(A,4,3)*ones(1,2);
    Ustar(2,:)=dete(A,3,2)/dete(A,3,4)*ones(2,1);
    Vstar(3,1)=-dete(A,[3,4],[2,3])*dete(A,4,1)/(dete(A,[3,4],[1,3])*dete(A,4,2));
    Vstar(4,1)=-dete(A,[3,4],[2,4])*dete(A,3,1)/(dete(A,[3,4],[1,4])*dete(A,3,2));
    Vstar(1,1)=dete(A,4,2)/dete(A,4,3);
    Vstar(2,2)=dete(A,3,2)/dete(A,3,4);
    vecXstar=Xstar(:);
    xstar=vecXstar(vecXstar~=0);
    %subs(solution.V4_1,a,A(:))-Vstar(4,1)
    norm(A-(Xstar*Ystar').*(Ustar*Vstar'))

    %% Things for the differential map: it is of full rank in nice cases
    % %grobner_basis=gbasis((system-a).',[x;u;y;v]);
    % nabla=jacobian(system-a,[x;y;u;v]);
    % rng(11)
    % % x=randn(length(x),1);
    % % y=randn(length(y),1);
    % % u=randn(length(u),1);
    % % v=randn(length(v),1);
    % % a=randn(length(a),1);
    % lengths=length(x)+length(u)+length(y)+length(v);
    % nablaeval=subs(nabla,[x;y;u;v],randn(lengths,1));
    % svd(nablaeval)
    % For this choice U(1:2,:)=[1,0;0,1];    V(1:2,:)=[1,1;1,1];
    % det(nabla(:,[1,3,4,6,7,9,10,12,13,15,16,18,19,20,23,24])) is not 0
    % Link for some theory
    %https://math.stackexchange.com/questions/1923855/
    % implicit-definition-of-the-intersection-of-n-hypersurfaces

    

function C=sym_fs(A,B)

    [n,p]=size(A);
    [~,q]=size(B);
    C=sym("C",[n,p*q]);
    for j=1:n
        C(j,:)=kron(A(j,:),B(j,:));
        %C(j,:)=reshape(A(j,:).'*B(j,:),p*q,1);
    end

end

%% Determinant of the submatrix of A obtained by removing rows and cols.
function a=dete(A,rows,cols)
    kr = true(size(A,1),1);
    kr(rows) = false;    % some rows to delete
    kc = true(1,size(A,2));
    kc(cols) = false;  % some columns to delete
    A=A(kr,kc);
    a=det(A);
end