function [Y1,Y2]=Upd_manifold(Y1,Y2,G,h)

    [n,r]=size(Y1);
    
    for i=1:n
        % Gi=reshape(G(i,:),r,r);
        % Gixi=Gi*Y1(i,:)';
        % Y1(i,:)=Y1(i,:)+h*(-Gi'*Y2(i,:)'+(Y2(i,:)*Gixi)*Y1(i,:)')';
        % Y1(i,:)=Y1(i,:)/norm(Y1(i,:));
        % Y2(i,:)=Y2(i,:)+h*(-Gixi)';

        Gi=reshape(G(i,:),r,r);
        Gixi=Gi*Y1(i,:)';
        uiGixi=Y2(i,:)*Gixi;
        Y1(i,:)=Y1(i,:)+h*(-Gi'*Y2(i,:)'+0.5*(uiGixi)*Y1(i,:)')';
        Y2(i,:)=Y2(i,:)+h*(-Gixi+0.5*(uiGixi)*Y2(i,:)')';
    end

end