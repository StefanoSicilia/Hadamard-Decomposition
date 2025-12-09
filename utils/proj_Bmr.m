function [A1,A2]=proj_Bmr(A)
%% proj_Bmr: 
% Projects the m-by-r^2 matrix A onto Bmr, that is the manifold of m-by-r^2
% matrices where each row is the vectorization of a rank-1 r-by-r matrix.
% It returns the factors such that proj_Bmr(A)=face_split(A1,A2), where A1
% and A2 are two m-by-r matrices.
    
    [m,rsq]=size(A);
    r=sqrt(rsq);
    A1=zeros(m,r);
    A2=zeros(m,r);
    for i=1:m
        [P,Sigma,Q]=svd(reshape(A(i,:),r,r));
        A1(i,:)=sqrt(Sigma(1,1))*Q(:,1)';
        A2(i,:)=sqrt(Sigma(1,1))*P(:,1)';
    end

end