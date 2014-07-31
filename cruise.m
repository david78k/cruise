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

% function load_data()
%------Load data of Challenge A---------%
handles.data = load('data_gui_08.mat');

handles.filter_a = handles.data.filter(1);
handles.filter_tau = handles.data.filter(2);
handles.gain_a = handles.data.gain(1);
handles.gain_b = handles.data.gain(2);
handles.gain_c = handles.data.gain(3);
handles.nu_0_threshold = handles.data.threshold(1);
handles.tau_threshold = handles.data.threshold(2);
handles.nu_1_threshold = handles.data.threshold(3);

handles.mean_g_e = [];
handles.std_g_e = [];
handles.mean_g_i = [];
handles.std_g_i = [];
handles.int_rel = [];
handles.coincidence = [];
handles.result = [];

handles.g_e = [];
handles.g_i = [];
handles.V = [];
handles.spiketimes = [];

handles.V_srm = [];
handles.filter = [];
handles.gain = [];
handles.threshold = [];
handles.coinc_spike = [];
%---------------------------------------%

