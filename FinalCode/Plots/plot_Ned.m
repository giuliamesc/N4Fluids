% Loading numerical solution
addpath('../Unstructured Mesh/');
eig_Ned = load('lambda_unstructured_Nedelec.dat')';

% Generating the first 40 exact eigenvalues
N = 40;
C=reshape((n2+m2.')',[],1)';
omega2 = sort(C);
omega2 = omega2(2:N+1);
x = [1:N];

% Plot
figure
plot(x,omega2,'bo','MarkerSize',10)
hold on
grid on
plot(x,eig_Ned(1:40),'xg','MarkerSize',10)
legend('Exact','Computed')
title('Nedelec, Unstructured Mesh')