lambdas = load('lambda1.dat')';
n2 = [0:7].^2;
m2 = [0:7].^2;

C=reshape((n2+m2.')',[],1)';
omega2 = sort(C);
omega2 = omega2(1:40);