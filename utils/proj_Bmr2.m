function [A1,A2]=proj_Bmr2(A)
%% proj_Bmr: orthogonal projection onto Bmr 
% Projects the m-by-r^2 matrix A onto Bmr, that is the manifold of m-by-r^2
% matrices where each row is the vectorization of a rank-1 r-by-r matrix.
% It returns the m-by-r factors A1 and A2 such that 
% norm(proj_Bmr(A)-face_split(A1,A2),'fro') is minimized.
    
    [m,rsq]=size(A);
    r=sqrt(rsq);

    % All-at-once rank-1 svd of the rows
    [P,Sigma,Q]=pagesvd(permute(reshape(A.',r,r,m),[1 2 3]));
    s=sqrt(squeeze(Sigma(1,1,:)));
    A1=s.*squeeze(Q(:,1,:)).';
    A2=s.*squeeze(P(:,1,:)).';

end