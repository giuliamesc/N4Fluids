n2 = [0:7].^2;
m2 = [0:7].^2;

N = 50;

C=reshape((n2+m2.')',[],1)';
omega2 = sort(C);
omega2 = omega2(2:N+1);
x = [1:N];
addpath('convtest/')
ref3 = load('convtestun1.dat')';
ref4 = load('convtestun2.dat')';
ref5 = load('convtestun3.dat')';
ref6 = load('convtestun4.dat')';
ref7 = load('convtestun5.dat')';
ref8 = load('convtestun6.dat')';
ref9 = load('convtestun7.dat')';

err = zeros(7,1);
err(1) = max(abs(omega2-ref3));
err(2) = max(abs(omega2-ref4));
err(3) = max(abs(omega2-ref5));
err(4) = max(abs(omega2-ref6));
err(5) = max(abs(omega2-ref7));
err(6) = max(abs(omega2-ref8));
err(7) = max(abs(omega2-ref9));

plot([3:9],err,'r*-')
hold on
grid on
yline(2*1e-2,'g')
legend('Error trend','Tolerance')

