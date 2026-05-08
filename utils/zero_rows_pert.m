function [W1,H1,W2,H2]=zero_rows_pert(W1,H1,W2,H2,theta)
%% zero_rows_pert: Perturbation of zero rows
% Replaces the zero rows of the matrices W1,H1,W2 and H2 with a random
% vector multiplied by the parameter theta.

    r1=size(W1,2);
    r2=size(W2,2);
    rng(1)
    vW1=vecnorm(W1,2,2)==0; W1(vW1,:)=theta*rand(sum(vW1),r1);
    vH1=vecnorm(H1,2,2)==0; H1(vH1,:)=theta*rand(sum(vH1),r1);
    vW2=vecnorm(W2,2,2)==0; W2(vW2,:)=theta*rand(sum(vW2),r2);
    vH2=vecnorm(H2,2,2)==0; H2(vH2,:)=theta*rand(sum(vH2),r2);

end