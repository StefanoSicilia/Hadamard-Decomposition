function x=hadLS(s,A,b)

    H=A'*diag(s.^2)*A;
    d=A'*(s.*b);
    %x=H\d;
    x=lsqminnorm(H,d);

end