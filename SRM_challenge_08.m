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

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SRM_challenge_08_OpeningFcn, ...
                   'gui_OutputFcn',  @SRM_challenge_08_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% Author: Skander Mensi, LCN, epfl. (skander.mensi@gmail.com)


% --- Executes just before SRM_challenge_08 is made visible.
function SRM_challenge_08_OpeningFcn(hObject, eventdata, handles, varargin)

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

%------Initialize Edit_text zone--------%
set(handles.mean_g_exc_edit,'String','----');
set(handles.std_g_exc_edit,'String','----');
set(handles.mean_g_inh_edit,'String','----');
set(handles.std_g_inh_edit,'String','----');
set(handles.coincidence_edit,'String','----');
set(handles.intrinsic_reliability_edit,'String','----');
set(handles.result_edit,'String','----');
set(handles.a_filter_edit,'String',num2str(handles.filter_a));
set(handles.tau_filter_edit,'String',num2str(handles.filter_tau/10));
set(handles.a_gain_edit,'String',num2str(handles.gain_a));
set(handles.b_gain_edit,'String',num2str(handles.gain_b));
set(handles.c_gain_edit,'String',num2str(handles.gain_c));
set(handles.nu_0_edit,'String',num2str(handles.nu_0_threshold));
set(handles.tau_edit,'String',num2str(handles.tau_threshold/10));
set(handles.nu_1_edit,'String',num2str(handles.nu_1_threshold));
%---------------------------------------%

%--------Initialize Slider--------------%
set(handles.a_filter_slider,'Value',handles.filter_a);
set(handles.tau_filter_slider,'Value',handles.filter_tau/10);
set(handles.a_gain_slider,'Value',handles.gain_a);
set(handles.b_gain_slider,'Value',handles.gain_b);
set(handles.c_gain_slider,'Value',handles.gain_c);
set(handles.nu_0_slider,'Value',handles.nu_0_threshold);
set(handles.tau_slider,'Value',handles.tau_threshold/10);
set(handles.nu_1_slider,'Value',handles.nu_1_threshold);
%---------------------------------------%

%-------Initialize Axes-----------------%
handles.temp = 0:0.1:(length(handles.data.v)-1)*0.1;
axes(handles.srm_vs_data_axes);
xlabel('time [ms]');
ylabel('U [mV]');
axis tight;
%---------------------------------------%

handles.output = hObject;
guidata(hObject, handles);


function varargout = SRM_challenge_08_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------Core-Function------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function load_popupmenu_Callback(hObject, eventdata, handles)

val = get(hObject, 'Value')-1;

if (val==0 || val==7)
    axes(handles.srm_vs_data_axes);
    cla;
    set(handles.mean_g_exc_edit,'String','----');
    set(handles.std_g_exc_edit,'String','----');
    set(handles.mean_g_inh_edit,'String','----');
    set(handles.std_g_inh_edit,'String','----');
    set(handles.coincidence_edit,'String','----');
    set(handles.intrinsic_reliability_edit,'String','----');
    set(handles.result_edit,'String','----');
else
    if(val<7)
        val = val;
    elseif(val>7)
        val = val-1;
    end
    handles.g_e = handles.data.g_e(:,val);
    handles.g_i = handles.data.g_i(:,val);
    handles.V = handles.data.v(:,val);
    handles.mean_g_e = mean(handles.g_e);
    handles.std_g_e = std(handles.g_e);
    handles.mean_g_i = mean(handles.g_i);
    handles.std_g_i = std(handles.g_i);
    handles.int_rel = handles.data.int_rel(val);
    handles.spiketimes = Extract_spiketimes(handles.V,10000);
    set(handles.mean_g_exc_edit,'String',num2str(handles.mean_g_e));
    set(handles.std_g_exc_edit,'String',num2str(handles.std_g_e));
    set(handles.mean_g_inh_edit,'String',num2str(handles.mean_g_i));
    set(handles.std_g_inh_edit,'String',num2str(handles.std_g_i));
    set(handles.intrinsic_reliability_edit,'String',num2str(handles.int_rel));
    set(handles.coincidence_edit,'String','----');
    set(handles.result_edit,'String','----');
    
    handles.temp = 0:0.1:(length(handles.V)-1)*0.1;
    axes(handles.srm_vs_data_axes);
    cla;
    plot(handles.temp, handles.V);
    xlabel('time [ms]');
    ylabel('Voltage [mV]');
    axis tight;

end

guidata(hObject, handles);


function reset_button_Callback(hObject, eventdata, handles)

clear handles.data;
set(handles.load_popupmenu,'val',1);

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

%------Initialize Edit_text zone--------%
set(handles.mean_g_exc_edit,'String','----');
set(handles.std_g_exc_edit,'String','----');
set(handles.mean_g_inh_edit,'String','----');
set(handles.std_g_inh_edit,'String','----');
set(handles.coincidence_edit,'String','----');
set(handles.intrinsic_reliability_edit,'String','----');
set(handles.result_edit,'String','----');
set(handles.a_filter_edit,'String',num2str(handles.filter_a));
set(handles.tau_filter_edit,'String',num2str(handles.filter_tau/10));
set(handles.a_gain_edit,'String',num2str(handles.gain_a));
set(handles.b_gain_edit,'String',num2str(handles.gain_b));
set(handles.c_gain_edit,'String',num2str(handles.gain_c));
set(handles.nu_0_edit,'String',num2str(handles.nu_0_threshold));
set(handles.tau_edit,'String',num2str(handles.tau_threshold/10));
set(handles.nu_1_edit,'String',num2str(handles.nu_1_threshold));
%---------------------------------------%

%--------Initialize Slider--------------%
set(handles.a_filter_slider,'Value',handles.filter_a);
set(handles.tau_filter_slider,'Value',handles.filter_tau/10);
set(handles.a_gain_slider,'Value',handles.gain_a);
set(handles.b_gain_slider,'Value',handles.gain_b);
set(handles.c_gain_slider,'Value',handles.gain_c);
set(handles.nu_0_slider,'Value',handles.nu_0_threshold);
set(handles.tau_slider,'Value',handles.tau_threshold/10);
set(handles.nu_1_slider,'Value',handles.nu_1_threshold);
%---------------------------------------%

%-------Initialize Axes-----------------%
handles.temp = 0:0.1:(length(handles.data.v)-1)*0.1;
axes(handles.srm_vs_data_axes);
cla;
xlabel('time [ms]');
ylabel('U [mV]');
axis tight;
%---------------------------------------%

handles.output = hObject;
guidata(hObject, handles);


function launch_button_Callback(hObject, eventdata, handles)

if(get(handles.load_popupmenu, 'Value')-1 == 0 || get(handles.load_popupmenu, 'Value')-1 == 7)
    set(handles.intrinsic_reliability_edit,'String','NOTHING');
    set(handles.coincidence_edit,'String','TO');
    set(handles.result_edit,'String','PLOT !!');
else
    handles.filter = [handles.filter_a handles.filter_tau];
    handles.gain = [handles.gain_a handles.gain_b handles.gain_c];
    handles.threshold = [handles.nu_0_threshold handles.tau_threshold handles.nu_1_threshold];
    [handles.V_srm handles.V_srm_spiketimes handles.thres] = ...
        adapt_gain_dif(handles.g_e,handles.g_i,handles.filter,handles.gain,10000,handles.threshold);
    [handles.coincidence handles.coinc_spike] = GamCoincFac(handles.V_srm_spiketimes, handles.spiketimes, 10000);
    handles.result = handles.coincidence/handles.int_rel;
    set(handles.coincidence_edit,'String',num2str(handles.coincidence));
    set(handles.result_edit,'String',num2str(handles.result));
    
    handles.temp = 0:0.1:(length(handles.data.v)-1)*0.1;
    axes(handles.srm_vs_data_axes);
    cla;
    hold on;
    plot(handles.temp, handles.V);
    %plot(handles.temp, handles.V_srm,'r');
    
    X = handles.temp(1:900);
    Y = handles.V_srm(1:900);
    p = plot(X(1:900),Y(1:900),'-r','EraseMode','none');
    
    for t=901:900:length(handles.temp);
        X = handles.temp(t-900:t);
        Y = handles.V_srm(t-900:t);
        set(p,'XData',X,'YData',Y) 
        drawnow
    end
    plot(handles.temp, handles.V_srm,'r');
    plot(handles.temp, handles.thres,'k');
    xlabel('time [ms]');
    ylabel('U [mV]');
    axis tight;
    if(isempty(handles.coinc_spike))
    else
        plot(handles.coinc_spike/10, 50,'kv');
        plot(10,55,'.');
        axis tight;
    end
    
end

handles.output = hObject;
guidata(hObject, handles);


function export_trace_button_Callback(hObject, eventdata, handles)

if(isempty(handles.V_srm) || isempty(handles.V))
    set(handles.intrinsic_reliability_edit,'String','NOTHING');
    set(handles.coincidence_edit,'String','TO');
    set(handles.result_edit,'String','PLOT !!');
else
    figure(2);
    hold on;
    plot(handles.temp, handles.V);
    plot(handles.temp, handles.V_srm,'r');
    plot(handles.temp, handles.thres,'k');
    xlabel('time [ms]');
    ylabel('membrane potential [mV]');
    axis tight;
    if(isempty(handles.coinc_spike))
    else
        plot(handles.coinc_spike/10, 50,'kv');
        axis tight;
    end
end

handles.output = hObject;
guidata(hObject, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------Edit text & Slider-------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function a_filter_slider_Callback(hObject, eventdata, handles)

handles.filter_a = get(handles.a_filter_slider,'Value');
set(handles.a_filter_edit,'String',num2str(handles.filter_a));

guidata(hObject, handles);


function a_filter_edit_Callback(hObject, eventdata, handles)

handles.filter_a = str2double(get(handles.a_filter_edit,'String'));

if(handles.filter_a < get(handles.a_filter_slider,'Min'))
    set(handles.a_filter_slider,'Min',handles.filter_a)
elseif(handles.filter_a > get(handles.a_filter_slider,'Max'))
    set(handles.a_filter_slider,'Max',handles.filter_a)
end

set(handles.a_filter_slider,'Value',handles.filter_a);

guidata(hObject, handles);


function tau_filter_slider_Callback(hObject, eventdata, handles)

handles.filter_tau = get(handles.tau_filter_slider,'Value');
set(handles.tau_filter_edit,'String',num2str(handles.filter_tau));

handles.filter_tau = handles.filter_tau*10;
guidata(hObject, handles);


function tau_filter_edit_Callback(hObject, eventdata, handles)

handles.filter_tau = str2double(get(handles.tau_filter_edit,'String'))*10;

if(handles.filter_tau/10 < get(handles.tau_filter_slider,'Min'))
    set(handles.tau_filter_slider,'Min',handles.filter_tau/10)
elseif(handles.filter_tau/10 > get(handles.tau_filter_slider,'Max'))
    set(handles.tau_filter_slider,'Max',handles.filter_tau/10)
end

set(handles.tau_filter_slider,'Value',handles.filter_tau/10);

guidata(hObject, handles);


function a_gain_slider_Callback(hObject, eventdata, handles)

handles.gain_a = get(handles.a_gain_slider,'Value');
set(handles.a_gain_edit,'String',num2str(handles.gain_a));

guidata(hObject, handles);


function a_gain_edit_Callback(hObject, eventdata, handles)

handles.gain_a = str2double(get(handles.a_gain_edit,'String'));

if(handles.gain_a < get(handles.a_gain_slider,'Min'))
    set(handles.a_gain_slider,'Min',handles.gain_a)
elseif(handles.gain_a > get(handles.a_gain_slider,'Max'))
    set(handles.a_gain_slider,'Max',handles.gain_a)
end

set(handles.a_gain_slider,'Value',handles.gain_a);

guidata(hObject, handles);


function b_gain_slider_Callback(hObject, eventdata, handles)

handles.gain_b = get(handles.b_gain_slider,'Value');
set(handles.b_gain_edit,'String',num2str(handles.gain_b));

guidata(hObject, handles);


function b_gain_edit_Callback(hObject, eventdata, handles)

handles.gain_b = str2double(get(handles.b_gain_edit,'String'));

if(handles.gain_b < get(handles.b_gain_slider,'Min'))
    set(handles.b_gain_slider,'Min',handles.gain_b)
elseif(handles.gain_b > get(handles.b_gain_slider,'Max'))
    set(handles.b_gain_slider,'Max',handles.gain_b)
end

set(handles.b_gain_slider,'Value',handles.gain_b);

guidata(hObject, handles);


function c_gain_slider_Callback(hObject, eventdata, handles)

handles.gain_c = get(handles.c_gain_slider,'Value');
set(handles.c_gain_edit,'String',num2str(handles.gain_c));

guidata(hObject, handles);


function c_gain_edit_Callback(hObject, eventdata, handles)

handles.gain_c = str2double(get(handles.c_gain_edit,'String'));

if(handles.gain_c < get(handles.c_gain_slider,'Min'))
    set(handles.c_gain_slider,'Min',handles.gain_c)
elseif(handles.gain_c > get(handles.c_gain_slider,'Max'))
    set(handles.c_gain_slider,'Max',handles.gain_c)
end

set(handles.c_gain_slider,'Value',handles.gain_c);

guidata(hObject, handles);


function nu_0_slider_Callback(hObject, eventdata, handles)

handles.nu_0_threshold = get(handles.nu_0_slider,'Value');
set(handles.nu_0_edit,'String',num2str(handles.nu_0_threshold));

guidata(hObject, handles);


function nu_0_edit_Callback(hObject, eventdata, handles)

handles.nu_0_threshold = str2double(get(handles.nu_0_edit,'String'));

if(handles.nu_0_threshold < get(handles.nu_0_slider,'Min'))
    set(handles.nu_0_slider,'Min',handles.nu_0_threshold)
elseif(handles.nu_0_threshold > get(handles.nu_0_slider,'Max'))
    set(handles.nu_0_slider,'Max',handles.nu_0_threshold)
end

set(handles.nu_0_slider,'Value',handles.nu_0_threshold);

guidata(hObject, handles);


function nu_1_slider_Callback(hObject, eventdata, handles)

handles.nu_1_threshold = get(handles.nu_1_slider,'Value');
set(handles.nu_1_edit,'String',num2str(handles.nu_1_threshold));

guidata(hObject, handles);


function nu_1_edit_Callback(hObject, eventdata, handles)

handles.nu_1_threshold = str2double(get(handles.nu_1_edit,'String'));

if(handles.nu_1_threshold < get(handles.nu_1_slider,'Min'))
    set(handles.nu_1_slider,'Min',handles.nu_1_threshold)
elseif(handles.nu_1_threshold > get(handles.nu_1_slider,'Max'))
    set(handles.nu_1_slider,'Max',handles.nu_1_threshold)
end

set(handles.nu_1_slider,'Value',handles.nu_1_threshold);

guidata(hObject, handles);


function tau_slider_Callback(hObject, eventdata, handles)

handles.tau_threshold = get(handles.tau_slider,'Value');
set(handles.tau_edit,'String',num2str(handles.tau_threshold));

handles.tau_threshold = handles.tau_threshold*10;
guidata(hObject, handles);


function tau_edit_Callback(hObject, eventdata, handles)

handles.tau_threshold = str2double(get(handles.tau_edit,'String'))*10;

if(handles.tau_threshold/10 < get(handles.tau_slider,'Min'))
    set(handles.tau_slider,'Min',handles.tau_threshold/10)
elseif(handles.tau_threshold/10 > get(handles.tau_slider,'Max'))
    set(handles.tau_slider,'Max',handles.tau_threshold/10)
end

set(handles.tau_slider,'Value',handles.tau_threshold/10);

guidata(hObject, handles);


function intrinsic_reliability_edit_Callback(hObject, eventdata, handles)


function coincidence_edit_Callback(hObject, eventdata, handles)


function result_edit_Callback(hObject, eventdata, handles)


function mean_g_exc_edit_Callback(hObject, eventdata, handles)


function std_g_exc_edit_Callback(hObject, eventdata, handles)


function mean_g_inh_edit_Callback(hObject, eventdata, handles)


function std_g_inh_edit_Callback(hObject, eventdata, handles)


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------Create-Function----------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


function std_g_inh_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function a_filter_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function a_filter_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_filter_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function tau_filter_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function a_gain_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function a_gain_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b_gain_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function b_gain_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function c_gain_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function c_gain_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function d_gain_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function d_gain_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function intrinsic_reliability_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function coincidence_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function result_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nu_0_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function nu_0_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nu_1_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function nu_1_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function tau_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function load_popupmenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mean_g_exc_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function std_g_exc_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mean_g_inh_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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

spiketimes_t(:,1) = spiketimes_t(:,1)-1;    %détecte le zéros crossing, sinon détecter le t d'après

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