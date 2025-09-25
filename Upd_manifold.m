function [X,U]=Upd_manifold(X,U,G,h)

    [n,r]=size(X);
    %W=face_split(X,U);
    %W=W+h*G;
    for i=1:n
        Gi=reshape(G(i,:),r,r);
        Gixi=Gi*X(i,:)';
        X(i,:)=X(i,:)+h*(-Gi'*U(i,:)'+(U(i,:)*Gixi)*X(i,:)')';
        X(i,:)=X(i,:)/norm(X(i,:));
        U(i,:)=U(i,:)+h*(-Gixi)';

        % Gi=reshape(G(i,:),r,r);
        % Gixi=Gi*X(i,:)';
        % uiGixi=U(i,:)*Gixi;
        % X(i,:)=X(i,:)+h*(-Gi'*U(i,:)'+0.85*(uiGixi)*X(i,:)')';
        % U(i,:)=U(i,:)+h*(-Gixi+0.15*(uiGixi)*U(i,:)')';

        % Gi=reshape(G(i,:),r,r);
        % UiGi=U(i,:)*Gi;
        % X(i,:)=X(i,:)+h*(-UiGi+(UiGi*X(i,:)')*X(i,:));
        % X(i,:)=X(i,:)/norm(X(i,:));
        % U(i,:)=U(i,:)+h*(-X(i,:)*Gi');

        % Gi=reshape(G(i,:),r,r);
        % Wi=reshape(W(i,:),r,r);
        % Wi=Wi+h*(-U(i,:)'*U(i,:)*Gi-Gi*X(i,:)'*X(i,:)+U(i,:)'*U(i,:)*Gi*X(i,:)'*X(i,:));
        % [Pi,Si,Qi]=svd(Wi);
        % U(i,:)=Si(1,1)*Pi(:,1)';
        % X(i,:)=Qi(:,1)';
    end

end