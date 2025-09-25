% Testing the decomposition in the case 4-by-4 with rank 2

    n=4;
    m=n;
    r=2;
    nmax=10;
    rng(1)

    A=randi(nmax,n,m);
    [P,S,Q]=svd(A);
    P=P*sqrt(S);
    Q=Q*sqrt(S);

    % a=sym("a",[r,1]);
    % b=sym("b",[r,1]);
    % c=sym("c",[r,1]);
    % d=sym("d",[r,1]);
    % alpha=sym("alpha",[n,1]);
    % beta=sym("beta",[m,1]);
    % rhs1=(alpha(1)*a+beta(1)*c).'*(alpha(2)*b+beta(2)*c);
    % rhs2=(alpha(3)*a+beta(3)*c).'*(alpha(4)*b+beta(4)*c);
    % lhs=a.'*b+c.'*d;
    % g=gbasis(rhs1*rhs2-lhs,[alpha,beta]);
    % h=solve(rhs1*rhs2-lhs,[alpha,beta]);

    % rng(1)
    % x1=randn(r,1); x1=x1/norm(x1);
    % y1=randn(r,1); y1=y1/norm(y1); 
    % x2=randn(r,1); x2=x2/norm(x2);
    % y2=randn(r,1); y2=y2/norm(y2);
    % a=x1'*y1; b=x1'*y2; c=x2'*y1; d=x2'*y2; e=x1'*x2; f=y1'*y2;
    % res1=cos(acos(c)-acos(b)-acos(a)-acos(d));
    % res2=cos(acos(b)+acos(c)+acos(a)+acos(d));
    % res3=cos(acos(a)-acos(b)-acos(c)+acos(d));
    % res4=cos(sign(a-e)*acos(a)+sign(b-f)*acos(b)+sign(c-f)*acos(c)+sign(d-e)*acos(d));
    % res5=cos(sign(f)*acos(a)+sign(e)*acos(b)+sign(f)*acos(c)+sign(e)*acos(d));
    % res6=cos(abs(pi/2-acos(a))+abs(pi/2-acos(b))-abs(pi/2-acos(c))-abs(pi/2-acos(d)));
    % res7=cos(sign(x1(1)-y1(1))*acos(a)+sign(x1(1)-y2(1))*acos(b)+...
    %     sign(x2(1)-y1(1))*acos(c)+sign(x2(1)-y2(1))*acos(d));
    % res8=atan2(det([x1,y1]),a)+atan2(det([y1,x2]),b)+...
    %     atan2(det([x2,y2]),d)+atan2(det([y2,x1]),c);
    % res9=1-cos(sign(det([x1,y2]))*acos(b)+sign(det([x2,y1]))*acos(c)+...
    %     sign(det([y1,x1]))*acos(a)+sign(det([y2,x2]))*acos(d));
    % res0=1-cos(sign(det([x1,y1]))*acos(a)+sign(det([y1,x2]))*acos(b)+...
    %     sign(det([x2,y2]))*acos(d)+sign(det([y2,x1]))*acos(c));
    % res=(a*b-sqrt((1-a^2)*(1-b^2)))*(c*d-sqrt((1-c^2)*(1-d^2)))-...
    % (sqrt(1-a^2)*b+sqrt(1-b^2)*a)*(sqrt(1-c^2)*d+sqrt(1-d^2)*c);
    % 
    % close all
    % plot([0,x1(1)],[0,x1(2)],'r-o')
    % hold on
    % plot([0,x2(1)],[0,x2(2)],'b-o')
    % hold on
    % plot([0,y1(1)],[0,y1(2)],'m-o')
    % hold on
    % plot([0,y2(1)],[0,y2(2)],'g-o')
    % hold on
    % plot(0,0,'ko')
    % legend('x1','x2','y1','y2','Location','best')


    rng(1)
    x1=randn(r,1); x1=x1/norm(x1);
    y1=randn(r,1); y1=y1/norm(y1); 
    x2=randn(r,1); x2=x2/norm(x2);
    y2=randn(r,1); y2=y2/norm(y2);
    a=x1'*y1; b=x1'*y2; c=x2'*y1; d=x2'*y2; 
    alpha11=acos(a); alpha12=acos(b); alpha21=acos(c); alpha22=acos(d);
    beta11=sign(det([x1,y1])); gamma12=sign(det([y1,x2]));
    beta22=sign(det([x2,y2])); gamma21=sign(det([y2,x1]));
    res=1-cos(beta11*alpha11+gamma12*alpha12+beta22*alpha22+gamma21*alpha21);

    err=(a*b-beta11*gamma12*sqrt((1-a^2)*(1-b^2)))*...
        (c*d-beta22*gamma21*sqrt((1-c^2)*(1-d^2)))-...
        (beta11*sqrt(1-a^2)*b+gamma12*sqrt(1-b^2)*a)*...
        (beta22*sqrt(1-c^2)*d+gamma21*sqrt(1-d^2)*c)-1;
   