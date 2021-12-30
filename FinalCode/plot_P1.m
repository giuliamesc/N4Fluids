N = 40;
x = [1:N];
% Generating the first 40 non-null exact eigenvalues
n2 = [0:7].^2;
m2 = [0:7].^2;
C=reshape((n2+m2.')',[],1)';
omega2 = sort(C);
omega2 = omega2(2:N+1);
% Loading eigenvalues computed with FeNiCs
eig_P1 = [1.000017, 1.000017, 2.000067, 4.000268, 4.000268, 5.000418, 5.000418, 5.999699, 8.001071, 9.001355, 9.001355, 10.001673, 10.001673, 13.002828, 13.002828, 14.995409, 14.995409, 16.004283, 16.004283, 17.004835, 17.004835, 18.005421, 20.006693, 20.006693, 23.995181, 25.010456, 25.010456, 25.010458, 25.010458, 26.011310, 26.011310, 29.014071, 29.014071, 29.973232, 29.973232, 32.017134, 34.019342, 34.019342, 36.021681, 36.021681];
% Vector of spurious eigenvalues and their positions
spurious_y = [6,15,15,24,29.97,29.97];
spurious_x = [8,16,17,25,34,35];
% Generating the figure
figure
plot(x,omega2,'bo','MarkerSize',10)
hold on
grid on
plot(x,eig_P1,'xg','MarkerSize',10)
plot(spurious_x,spurious_y,'xr','MarkerSize',10)
legend('Exact','Computed','Spurious')
title('P1, Criss-cross Mesh')