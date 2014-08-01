function cruise()
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
% step(r*T,t);
% axis([0 20 0 10])

Kp = 1000;
Ki = 50;
Kd = 1;
C = pid(Kp,Ki,Kd)

T = feedback(C*P_cruise,1)
% step(r*T,t)
stepinfo(r*T)

% handles = []
% load_data (handles)

% handles

% function load_data(handles)
%------Load data of Challenge A---------%
handles.data = load('data_gui_08.mat')

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
handles.coinc_spike = []
%---------------------------------------%

% call SRM 
handles.filter = [handles.filter_a handles.filter_tau];
handles.gain = [handles.gain_a handles.gain_b handles.gain_c];
handles.threshold = [handles.nu_0_threshold handles.tau_threshold handles.nu_1_threshold];
[handles.V_srm handles.V_srm_spiketimes handles.thres] = ...
    adapt_gain_dif(handles.g_e,handles.g_i,handles.filter,handles.gain,10000,handles.threshold);
[handles.coincidence handles.coinc_spike] = GamCoincFac(handles.V_srm_spiketimes, handles.spiketimes, 10000);
% handles.result = handles.coincidence/handles.int_rel;
% set(handles.coincidence_edit,'String',num2str(handles.coincidence));
% set(handles.result_edit,'String',num2str(handles.result));

handles.temp = 0:0.1:(length(handles.data.v)-1)*0.1
% axes(handles.srm_vs_data_axes);
cla;
hold on;
% plot(handles.temp, handles.V);
%plot(handles.temp, handles.V_srm,'r');

X = handles.temp(1:900);
% Y = handles.V_srm(1:900);
% p = plot(X(1:900),Y(1:900),'-r','EraseMode','none');

for t=901:900:length(handles.temp);
    X = handles.temp(t-900:t);
%     Y = handles.V_srm(t-900:t);
%     set(p,'XData',X,'YData',Y) 
    drawnow
end
% plot(handles.temp, handles.V_srm,'r');
% plot(handles.temp, handles.thres,'k');
xlabel('time [ms]');
ylabel('U [mV]');
axis tight;
if(isempty(handles.coinc_spike))
else
    plot(handles.coinc_spike/10, 50,'kv');
    plot(10,55,'.');
    axis tight;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------MODEL SRM----------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [U_sub spiketimes thres] = adapt_gain_dif(g_e,g_i,filter,gain,sampling_freq,threshold)
% [handles.V_srm handles.V_srm_spiketimes handles.thres] =
%   adapt_gain_dif(handles.g_exc,handles.g_inh,handles.filter,handles.gain,10000,handles.threshold)
%   Detailed explanation goes here

a = -filter(1);
tau = filter(2);
E_rest = -63.4843;
E_e = -10;
E_i = -70;
nu_0 = threshold(1);
tau_1 = threshold(2);
nu_1 = threshold(3);

    U_sub = zeros(length(g_e),1);
    V = E_rest;
    U_sub(1) = V;
    thres = zeros(length(g_e),1);
    thres(1) = nu_0;

    for t = 1:length(g_e)
        I = g_e(t)*(V-E_e) + g_i(t)*(V-E_i);
        dV = -V/tau + a*I + E_rest/tau;
        V = V + dV;
        U_sub(t) = V;
    end

    U_sub = U_sub(1:length(g_e));

    U_sub = gain(1)*U_sub.^3 + gain(2)*U_sub.^2 + gain(3)*U_sub;

    temp_spike = 0;
    t_refr = 202;
    U_threshold = nu_0;
    for j = 1:length(U_sub)

        if(temp_spike(end) == j+9)
            dthreshold = ((U_threshold - nu_0)/-tau_1) + nu_1;
        else
            dthreshold = ((U_threshold - nu_0)/-tau_1);
        end

        U_threshold = U_threshold + dthreshold;
        thres(j) = U_threshold;

        if(U_sub(j) >= U_threshold && t_refr > 200)
            U_sub(j+10) = 35;
            t_refr = 0;
            temp_spike(end+1) = j+10;
        else
            t_refr = t_refr + 1;
        end
    end

    spike = Extract_spiketimes(U_sub, sampling_freq);
    spiketimes = spike(:,1);


function [spiketimes] = Extract_spiketimes(voltage, sampling_freq)
%   sampling_freq is used to caculate the period refractory around
%   0.8 ms to detect any spikes.
%   limit_t_refr = floor(0.8*(sampling_freq*1e-3))
%   extract the spiketimes a trace. spiketimes(:,1) = spiketimes
%   the spiketimes are detected in zéros upward crossing
%   spiketimes(:,2) = max voltage of the spike, used for the plotter the detection
%   of spikes.
%   spiketimes(:,3) = ISI, useful for eta dependent the last isi;

k=1;                                                %counts the number of spike
voltage_prime = [0;diff(voltage)];                  %derived from voltage
limite = 0;                                         %limit from which it takes a spike = 0
limit_t_refr = floor(0.8*(sampling_freq*1e-3));     %see above
t_refr = limit_t_refr + 1;                          %allows to detect a spike in t=1
spiketimes_t(1,1) = 0;
spiketimes_t(1,2) = 0;
spiketimes_t(1,3) = 0;

for i=1:length(voltage)-6          %loop the voltage
    if(voltage(i) >= limite && voltage_prime(i) > 0 && t_refr >= limit_t_refr)
        k = k+1;
        t_refr = 0;
        spiketimes_t(k,2) = max(voltage(i:i+5)); % 1 takes on the max [ms]
        spiketimes_t(k,1) = i;
        spiketimes_t(k,3) = spiketimes_t(k,1)-spiketimes_t(k-1,1);
    else
        t_refr = t_refr + 1;
    end
end

spiketimes_t(:,1) = spiketimes_t(:,1)-1;    %detect the zéros crossing, if the detection t after

spiketimes(:,1) = spiketimes_t(2:end,1);
spiketimes(:,2) = spiketimes_t(2:end,2);
spiketimes(:,3) = spiketimes_t(2:end,3);


%
%    G = GamCoincFac(PredSpkTrain,TargetSpkTrain, SamplingFreq)
%       Calculates the Gamma factor for the Target Spike Train
%       TargetSpkTrain and the Modeled Spike Train PredSpkTrain.
%       SamplingFreq the sampling frequency in Hz, and the spike trains are
%       a list of spike indices (spike times = SpkTrain / SamplingFreq).  
%   
%   See
%       Kistler et al, Neural Comp 9:1015-1045 (1997)
%       Jolivet et al, J Neurophysiol 92:959-976 (2004)
%   for further details
%
%
%           - Renaud Jolivet 2007.05.
%           - modified by R. Naud 2007.09. (siingularity and break in loop
%          and empty PredSpkTrain handling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [G coinc_spike] = GamCoincFac(PredSpkTrain,TargetSpkTrain,SamplingFreq)
% [handles.coincidence handles.coinc_spike] = GamCoincFac(handles.V_srm_spiketimes, handles.spiketimes, 10000);

coinc_spike = [];

if isempty(PredSpkTrain), G = 0; return, end


% Some parameters
DeltaWindow     =   2e-03;                      % [sec]         % half-width of coincidence detection (2 msec)
DeltaBins       =   DeltaWindow*SamplingFreq;   % [time bins]   % half-width of coincidence detection (2 msec)
NSpikesPred     =   length(PredSpkTrain);
NSpikesTarget   =   length(TargetSpkTrain);
%
% Compute frequencies, normalisation and average coincidences 
FreqPred        =   SamplingFreq*(NSpikesPred-1)/max((PredSpkTrain(NSpikesPred)-PredSpkTrain(1)),1);
NCoincAvg       =   2*DeltaWindow*NSpikesTarget*FreqPred;
NNorm           =   abs(1-2*FreqPred*DeltaWindow);
%
% Compute the gamma coincidence factor
NCoinc          =   0; 
i               =   1;
while i <= NSpikesTarget
    j=1;
    while j <= NSpikesPred
        if abs(PredSpkTrain(j)-TargetSpkTrain(i)) <= DeltaBins
            NCoinc  =   NCoinc+1;
            i       =   i+1;
            coinc_spike(end+1) = PredSpkTrain(j);
            if i> NSpikesTarget, break, end
        end
        j=j+1;
    end
    i=i+1;
end

G               =	(NCoinc-NCoincAvg)/(1/2*(NSpikesPred+NSpikesTarget))*1/NNorm;

if(G<0)
	G = 0;
end

