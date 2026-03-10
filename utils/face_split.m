function C=face_split(A,B)
%% face_split:
% Computes the face_split product of A and B. It assumes that A and B have
% the same number of rows. It avoids the loop and kron product of the
% mathematical definition, that is
% C=zeros(m,n*p);
% for j=1:m
%     C(j,:)=kron(A(j,:),B(j,:));
% end
% and it implements the product via Matlab reshapings.
    
    [m,r1]=size(A);
    [~,r2]=size(B);
    C=reshape(permute(reshape(A,m,r1,1).*reshape(B,m,1,r2),[1 3 2]),m,[]);

end