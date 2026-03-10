function [A1,A2]=projBmr_loop(A)
%% projBmr_loop: orthogonal projection onto Bmr 
% Projects the m-by-r^2 matrix A onto Bmr, that is the manifold of m-by-r^2
% matrices where each row is the vectorization of a rank-1 r-by-r matrix.
% It returns the m-by-r factors A1 and A2 such that 
% norm(projBmr(A)-face_split(A1,A2),'fro') is minimized.
% Same as projBmr_noloop, but it uses a loop.
    
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