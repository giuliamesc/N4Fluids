lambda1 = load('lambda1.dat')';
lambda2 = load('lambda2.dat');
lambda3 = load('lambda3.dat');
lambda4 = load('lambda4.dat');

n2 = [0:7].^2;
m2 = [0:7].^2;

N = 30;

C=reshape((n2+m2.')',[],1)';
omega2 = sort(C);
omega2 = omega2(2:N+1);
x = [1:N];

figure
plot(x,omega2,'bo','MarkerSize',10)
hold on
grid on
plot(x,lambda1,'xr','MarkerSize',10)
legend('exact','computed')
title('Nedelec, Triangular Mesh')

figure
plot(x,omega2,'bo','MarkerSize',10)
hold on
grid on
plot(x,lambda2,'xr','MarkerSize',10)
legend('exact','computed')
title('P1, Triangular Mesh')

figure
plot(x,omega2,'bo','MarkerSize',10)
hold on
grid on
plot(x,lambda3,'xr','MarkerSize',10)
legend('exact','computed')
title('Nedelec, Unstructured Mesh')

figure
plot(x,omega2,'bo','MarkerSize',10)
hold on
grid on
plot(x,lambda4,'xr','MarkerSize',10)
legend('exact','computed')
title('P1, Unstructured Mesh')