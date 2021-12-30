n2 = [0:6].^2;
m2 = [0:6].^2;

N = 48;

nbRefs = 6

C=reshape((n2+m2.')',[],1)';
omega2 = sort(C);
omega2 = omega2(2:N+1);
x = [1:N];
ref3 = load('lambda_unstructured_Nedelec_convergence_3.dat')';
ref4 = load('lambda_unstructured_Nedelec_convergence_4.dat')';
ref5 = load('lambda_unstructured_Nedelec_convergence_5.dat')';
ref6 = load('lambda_unstructured_Nedelec_convergence_6.dat')';
ref7 = load('lambda_unstructured_Nedelec_convergence_7.dat')';
ref8 = load('lambda_unstructured_Nedelec_convergence_8.dat')';

ref3 = sort(ref3);
ref4 = sort(ref4);
ref5 = sort(ref5);
ref6 = sort(ref6);
ref7 = sort(ref7);
ref8 = sort(ref8);

err = zeros(nbRefs,1);
err(1) = max(abs(omega2-ref3));
err(2) = max(abs(omega2-ref4));
err(3) = max(abs(omega2-ref5));
err(4) = max(abs(omega2-ref6));
err(5) = max(abs(omega2-ref7));
err(6) = max(abs(omega2-ref8));

plot([3:2 + nbRefs],err,'r*-')
hold on
grid on
legend('Error trend','Tolerance')

