function varargout = SRM_challenge_08(varargin)
% SRM_CHALLENGE_08 M-file for SRM_challenge_08.fig
%
%      to run this programm: >SRM_challenge_08
%
%      SRM_CHALLENGE_08, by itself, creates a new SRM_CHALLENGE_08 or raises the existing
%      singleton*.
%
%      H = SRM_CHALLENGE_08 returns the handle to a new SRM_CHALLENGE_08 or the handle to
%      the existing singleton*.
%
%      SRM_CHALLENGE_08('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SRM_CHALLENGE_08.M with the given input arguments.
%
%      SRM_CHALLENGE_08('Property','Value',...) creates a new SRM_CHALLENGE_08 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SRM_challenge_08_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SRM_challenge_08_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 04-Mar-2008 12:31:28

% [handles.V_srm handles.V_srm_spiketimes handles.thres] = ...
%         adapt_gain_dif(handles.g_e,handles.g_i,handles.filter,handles.gain,10000,handles.threshold);
    
% [handles.coincidence handles.coinc_spike] = ...
% GamCoincFac(handles.V_srm_spiketimes, handles.spiketimes, 10000);

% Begin initialization code - DO NOT EDIT
% gui_Singleton = 1;
% gui_State = struct('gui_Name',       mfilename, ...
%                    'gui_Singleton',  gui_Singleton, ...
%                    'gui_OpeningFcn', @SRM_challenge_08_OpeningFcn, ...
%                    'gui_OutputFcn',  @SRM_challenge_08_OutputFcn, ...
%                    'gui_LayoutFcn',  [] , ...
%                    'gui_Callback',   []);
% if nargin && ischar(varargin{1})
%     gui_State.gui_Callback = str2func(varargin{1});
% end
% 
% if nargout
%     [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
% else
%     gui_mainfcn(gui_State, varargin{:});
% end
% End initialization code - DO NOT EDIT

% Author: Skander Mensi, LCN, epfl. (skander.mensi@gmail.com)


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------MODEL SRM----------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [U_sub spiketimes thres] = adapt_gain_dif(g_e,g_i,filter,gain,sampling_freq,threshold)
%   [spiketimes] = adapt_gain_dif(g_exc,g_inh,filter,gain,10000,threshold)
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
%   sampling_freq est utile pour caculer la période réfractaire, environ
%   0.8 ms pour détecter tout les spikes.
%   limit_t_refr = floor(0.8*(sampling_freq*1e-3))
%   extrait les spiketimes d'un trace. spiketimes(:,1) = spiketimes
%   les spiketimes sont detecté en zéros upward crossing
%   spiketimes(:,2) = max voltage du spike, utile pour plotter la détection
%   des spikes.
%   spiketimes(:,3) = ISI, utile pour eta dépendant du last isi;

k=1;                                                %compte le nombre de spike
voltage_prime = [0;diff(voltage)];                  %dérivée de voltage
limite = 0;                                         %limite à partir de laquelle on prend un spike = 0
limit_t_refr = floor(0.8*(sampling_freq*1e-3));     %voir au dessus
t_refr = limit_t_refr + 1;                          %permet de detecter un spike en t=1
spiketimes_t(1,1) = 0;
spiketimes_t(1,2) = 0;
spiketimes_t(1,3) = 0;

for i=1:length(voltage)-6          %parcours le voltage
    if(voltage(i) >= limite && voltage_prime(i) > 0 && t_refr >= limit_t_refr)
        k = k+1;
        t_refr = 0;
        spiketimes_t(k,2) = max(voltage(i:i+5)); % prend le max sur 1 [ms]
        spiketimes_t(k,1) = i;
        spiketimes_t(k,3) = spiketimes_t(k,1)-spiketimes_t(k-1,1);
    else
        t_refr = t_refr + 1;
    end
end

spiketimes_t(:,1) = spiketimes_t(:,1)-1;    %détecte le zéros crossing, sinon détecter le t d'apr?s

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
else
	G = G;
end