function x=hadLS(s,A,b)
%% hadLS:
% Least square solution required for the BCD method. For the details see  
% [WVG25] S. Wertz, A. Vandaele, N.Gillis, 
% Efficient algorithms for the Hadamard decomposition, 2025.

    H=A'*diag(s.^2)*A;
    d=A'*(s.*b);
    x=lsqminnorm(H,d); %x=H\d;

end