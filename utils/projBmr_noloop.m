function [A1,A2]=projBmr_noloop(A,r1,r2)
%% projBmr_noloop: orthogonal projection onto Bmr 
% Projects the m-by-r^2 matrix A onto Bmr, that is the manifold of m-by-r^2
% matrices where each row is the vectorization of a rank-1 r-by-r matrix.
% It returns the m-by-r factors A1 and A2 such that 
% norm(proj_Bmr(A)-face_split(A1,A2),'fro') is minimized.
% Same as projBmr_loop, but it uses Matlab reshapings to avoid the loop.
    
    m=size(A,1);

    % All-at-once rank-1 svd of the rows
    [P,Sigma,Q]=pagesvd(permute(reshape(A.',r2,r1,m),[1 2 3]));
    s=sqrt(squeeze(Sigma(1,1,:)));
    A1=s.*squeeze(Q(:,1,:)).';
    A2=s.*squeeze(P(:,1,:)).';

end