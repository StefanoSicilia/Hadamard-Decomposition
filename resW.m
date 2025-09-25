function a=resW(R,A,H)

% R={{x1,...,xn},{u1,...,un}}

    [n,~]=size(A);
    [~,rsq]=size(H);
    W=zeros(n,rsq);
    for i=1:n
        W(i,:)=reshape(R{2}{i}'*R{1}{i},rsq,1)';
    end
    a=norm(A-W*H','fro')^2;
    
end