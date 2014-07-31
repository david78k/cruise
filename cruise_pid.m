% rise time < 5s
% overshoot < 10%
% steady state error < 2%

% vehicle mass in kg
m = 1000; 
% damping coefficient in N.s/m
b = 50; 
% reference speed in m/s
r = 10; 

s = tf('s');
P_cruise = 1/(m*s + b);

% PID control
Kp = 1;
Ki = 1;
Kd = 1;
C = pid(Kp,Ki,Kd);

T = feedback(C*P_cruise,1);
t = 0:0.1:20;
step(r*T,t);
% axis([0 20 0 10])

Kp = 1000;
Ki = 50;
Kd = 1;
C = pid(Kp,Ki,Kd)

T = feedback(C*P_cruise,1)
step(r*T,t)
stepinfo(r*T)

% Proportional control
% Kp = 1;
% Ki = 1;
% Kd = 1;
% 
% s = tf('s');
% C = Kp + Ki/s + Kd*s
% C = pid(Kp,Ki,Kd)
% 
% Kp = 100;
% C = pid(Kp);
% 
% T = feedback(C*P_cruise,1)
% 
% t = 0:0.1:20;
% step(r*T,t)
% axis([0 20 0 10])
% 
% % 
% Kp = 5000;
% C = pid(Kp);
% T = feedback(C*P_cruise,1);
% 
% step(r*T,t)
% axis([0 20 0 10])
% 
% % PI control
% Kp = 600;
% Ki = 1;
% C = pid(Kp,Ki);
% 
% T = feedback(C*P_cruise,1);
% 
% step(r*T,t)
% axis([0 20 0 10])
% 
% %
% Kp = 800;
% Ki = 40;
% C = pid(Kp,Ki);
% 
% T = feedback(C*P_cruise,1);
% 
% step(r*T,t)
% axis([0 20 0 10])

