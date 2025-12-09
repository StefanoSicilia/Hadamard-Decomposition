function C=face_split(A,B)
    
    [n,p]=size(A);
    [m,q]=size(B);
    % if m~=n
    %     error('Incompatible sizes of the matrices.')
    % end
    C=zeros(n,p*q);
    for j=1:n
        C(j,:)=kron(A(j,:),B(j,:));
    end
end