u = [Open High Low];
t = Date;
zero_matrix = ones(length(t),1);
y = [Close];
data = iddata(y,u,1);
opt = ssestOptions('InitialState','auto','Focus','simulation','Display','on');
[sys,Parameters] = ssest(data(1:2333),3,'Ts',1,'Form','free','DisturbanceModel','none',opt);
A=sys.A;
B=sys.B;
C=sys.C;

Q=eye(3);
R=1;
n = length(t);
P = B*Q*B';         % Initial error covariance
x = zeros(3,1);     % Initial condition on the state
%ye = zeros(length(t),3);
%ycov = zeros(length(t),3);
%y1 = zeros(length(t),3);

for i = 1:length(t)
  % Measurement update
  Mn = P*C'/(C*P*C'+R);
  x = x + Mn*(y(i)-C*x);   % x[n|n]
  P = (eye(3)-Mn*C)*P; 
  % P[n|n]
  ye(i) =  C*x;
  % errcov(i) = C*P*C';
  % Time update
  H=[u(i);u(i+1);u(i+2)];
  x = A*x + B*H;
  y1(i)=C*x;% x[n+1|n]
  P = A*P*A' + B*Q*B';     % P[n+1|n]
end


for i=1:length(t)
    yr(i,1)=ye(1,i);
    yp(i,1)=y1(1,i);
end

subplot(211), plot(t,y,'--',t,ye,'-')
title('Time-varying Kalman filter response')
xlabel('time'), ylabel('Output')
subplot(212), plot(t,yp,'-')
title('Time-varying Kalman filter response')
xlabel('time'), ylabel('Output')

%subplot(212), plot(t,y-yv,'-.',t,y-ye,'-')
%xlabel('No. of samples'), ylabel('Output')

%subplot(211)
%plot(t,errcov), ylabel('Error covar')
EstErr = y - yr;
EstErrCov = sum(EstErr.*EstErr)/length(EstErr)

y1(n)
