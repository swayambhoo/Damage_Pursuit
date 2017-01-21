function varargout = Damage_Pursuit_Gui(varargin)
% Damage_Pursuit_Gui MATLAB code for Damage_Pursuit_Gui.fig
%      Damage_Pursuit_Gui, by itself, creates a new Damage_Pursuit_Gui or raises the existing
%      singleton*.
%
%      H = Damage_Pursuit_Gui returns the handle to a new Damage_Pursuit_Gui or the handle to
%      the existing singleton*.
%
%      Damage_Pursuit_Gui('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Damage_Pursuit_Gui.M with the given input arguments.
%
%      Damage_Pursuit_Gui('Property','Value',...) creates a new Damage_Pursuit_Gui or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Damage_Pursuit_Gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Damage_Pursuit_Gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Damage_Pursuit_Gui

% Last Modified by GUIDE v2.5 12-Apr-2016 11:26:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Damage_Pursuit_Gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Damage_Pursuit_Gui_OutputFcn, ...
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


% --- Executes just before Damage_Pursuit_Gui is made visible.
function Damage_Pursuit_Gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Damage_Pursuit_Gui (see VARARGIN)



% Choose default command line output for Damage_Pursuit_Gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes Damage_Pursuit_Gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Damage_Pursuit_Gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function tau1Val_Callback(hObject, eventdata, handles)
% hObject    handle to tau1Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau1Val as text
%        str2double(get(hObject,'String')) returns contents of tau1Val as a double



% --- Executes during object creation, after setting all properties.
function tau1Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau1Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = waitbar(0,'Loading data and configurations ...');

filename = strcat(handles.pathname,handles.filename);
Data = load(filename);
FldNm = fieldnames(Data);
handles.Data = getfield(Data,FldNm{1});

Spt_Rsl = get(handles.SptRsl,'string');
Spt_Rsl = str2double(Spt_Rsl);

T = size(handles.Data,3);
M = size(handles.Data,1);
N = size(handles.Data,2);

if (~isnan(Spt_Rsl)) 
    if (rem(M,Spt_Rsl)> 0) || (rem(N,Spt_Rsl)> 0)
        
        warning('The frames of the data are cropped to make the dimensions integer multiples of Spt_Rsl.');
        M = floor(M/Spt_Rsl) * Spt_Rsl;
        N = floor(N/Spt_Rsl) * Spt_Rsl;
        
        handles.Data = handles.Data(1:M,1:N,:);
    end
end

Lx = get(handles.LxVal,'string');
Lx = str2double(Lx);

if isnan(Lx) || (Lx <= 0)
    error('Please enter a valid value for Lx!');
end

Ly = get(handles.LyVal,'string');
Ly = str2double(Ly);

if isnan(Ly) || (Ly <= 0)
    error('Please enter a valid value for Ly!');
end

tau1 = get(handles.tau1Val,'string');
tau1 = str2double(tau1);

if isnan(tau1) || (tau1 <= 0)
    error('Please enter a valid value for tau1!');
end

tau2 = get(handles.tau2Val,'string');
tau2 = str2double(tau2);

if isnan(tau2) || (tau2 <= 0)
    error('Please enter a valid value for tau2!');
end

smtDictSel = get(handles.smtDictSel,'SelectedObject');
smtDictSel = get(smtDictSel,'string');
D1_Type = strcmp(smtDictSel,'DFT');

sprsDictSel = get(handles.sprsDictSel,'SelectedObject');
sprsDictSel = get(sprsDictSel,'string');
D2_Type = strcmp(sprsDictSel,'Marr Wavelet');

Tmax = 1000;

[startFrame, endFrame] = selectFrames(handles);


Y = handles.Data(:,:,startFrame:endFrame);
waitbar(0,h,'Running the algorithm ...');

if D2_Type
    
    sigma = get(handles.sigmaVal,'string');
    sigma = str2double(sigma);
    
    if isnan(sigma) || (sigma <= 0)
       error('Please enter a valid value for sigma!'); 
    end
    
    thr = 1e-5;
    
    for i = 1 : length(sigma)
        
        [d,mask] = MarrWvlt_Dct(Lx,Ly,M,N,sigma(i),thr(i));
        dMtx(:,i) = d;
        maskMtx(:,i) = mask;
        
    end
    
    handles.dMtx = dMtx;
    handles.maskMtx = maskMtx;
    
    
    [X1hat , X2hat, Obj , time ] = Damage_Pursuit...
    (Y,'sgm',sigma,'thr',thr,'Lx',Lx,'Ly',Ly,...
     'D1_Type',D1_Type,'D2_Type',D2_Type,'tau1',tau1,'tau2',tau2);

else
    
    Spt_Rsl = get(handles.SptRsl,'string');
    Spt_Rsl = str2double(Spt_Rsl);
    
    if isnan(Spt_Rsl)
        Spt_Rsl = 1;
    end
    
    if ~(Spt_Rsl <= min(M,N))
        error('Please enter a valid value for spatial resolution!')
    end
    
    [X1hat , X2hat, Obj , time ] = Damage_Pursuit...
    (Y,'Lx',Lx,'Ly',Ly,...
     'D1_Type',D1_Type,'D2_Type',D2_Type,'tau1',tau1,'tau2',...
     tau2,'Max_Itr',Tmax, 'Spt_Rsl',Spt_Rsl);


end

handles.output = [X1hat;X2hat];
handles.outputIndctr = [startFrame:endFrame];
% Update handles structure
guidata(hObject, handles);

waitbar(1,h,'Successful Implementation');
close(h);

set(handles.startPlot,'String',num2str(startFrame));
updateAxes(handles,startFrame);


function tau2Val_Callback(hObject, eventdata, handles)
% hObject    handle to tau2Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau2Val as text
%        str2double(get(hObject,'String')) returns contents of tau2Val as a double


% --- Executes during object creation, after setting all properties.
function tau2Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau2Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadFileButton.
function loadFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.filename, handles.pathname] = uigetfile({'*.mat'},'File Selector');
guidata(hObject, handles);



function LyVal_Callback(hObject, eventdata, handles)
% hObject    handle to LyVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LyVal as text
%        str2double(get(hObject,'String')) returns contents of LyVal as a double


% --- Executes during object creation, after setting all properties.
function LyVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LyVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LxVal_Callback(hObject, eventdata, handles)
% hObject    handle to LxVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LxVal as text
%        str2double(get(hObject,'String')) returns contents of LxVal as a double


% --- Executes during object creation, after setting all properties.
function LxVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LxVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigmaVal_Callback(hObject, eventdata, handles)
% hObject    handle to sigmaVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmaVal as text
%        str2double(get(hObject,'String')) returns contents of sigmaVal as a double


% --- Executes during object creation, after setting all properties.
function sigmaVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes during object creation, after setting all properties.
function radiobutton3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function smtDictSel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smtDictSel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function startFrame_Callback(hObject, eventdata, handles)
% hObject    handle to startFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFrame as text
%        str2double(get(hObject,'String')) returns contents of startFrame as a double


% --- Executes during object creation, after setting all properties.
function startFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endFrame_Callback(hObject, eventdata, handles)
% hObject    handle to endFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrame as text
%        str2double(get(hObject,'String')) returns contents of endFrame as a double


% --- Executes during object creation, after setting all properties.
function endFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function numFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function startPlot_Callback(hObject, eventdata, handles)
% hObject    handle to startPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

startPlot = get(handles.startPlot,'string');
startPlot = str2double(startPlot);


updateAxes(handles,startPlot);
% Hints: get(hObject,'String') returns contents of startPlot as text
%        str2double(get(hObject,'String')) returns contents of startPlot as a double


% --- Executes during object creation, after setting all properties.
function startPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotNext.
function plotNext_Callback(hObject, eventdata, handles)
% hObject    handle to plotNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[startFrame, endFrame] = selectFrames(handles);

startPlot = get(handles.startPlot,'string');
startPlot = str2double(startPlot);

if startPlot < endFrame
    startPlot = startPlot + 1;
end

updateAxes(handles,startPlot);
set(handles.startPlot, 'String', num2str(startPlot));


% --- Executes on button press in plotBack.
function plotBack_Callback(hObject, eventdata, handles)
% hObject    handle to plotBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[startFrame, endFrame] = selectFrames(handles);

startPlot = get(handles.startPlot,'string');
startPlot = str2double(startPlot);

if startPlot > startFrame
    startPlot = startPlot - 1;
end

updateAxes(handles,startPlot);
set(handles.startPlot, 'String', num2str(startPlot));


function [startFrame, endFrame] = selectFrames(handles)

numFrames = get(handles.numFrames,'SelectedObject');
numFrames = get(numFrames,'string');

if strcmp(numFrames,'All Frames')
    
    startFrame = 1;
    endFrame = size(handles.Data,3);

else
    
    startFrame = get(handles.startFrame,'string');
    startFrame = str2double(startFrame);
    
    if isnan(startFrame) || ~(1 <= startFrame) || ~(startFrame <= size(handles.Data,3))
        error('Please enter a valid value for the starting frame!');
    end
    
    endFrame = get(handles.endFrame,'string');
    endFrame = str2double(endFrame);
    
    if isnan(endFrame) || ~(1 <= endFrame) || ~(endFrame <= size(handles.Data,3))
        error('Please enter a valid value for the ending frame!');
    end
    
    if startFrame > endFrame
       error('The starting frame should have smaller frame number.');  
    end
    
end


function updateAxes(handles,startPlot)

[startFrame, endFrame] = selectFrames(handles);

if (startPlot < startFrame) || (startPlot > endFrame)
    error('Please enter a valid frame number to plot!');
end

Z = handles.Data(:,:,startPlot);

M = size(handles.Data,1);
N = size(handles.Data,2);

Xhat = handles.output;
indctr = find(handles.outputIndctr == startPlot);
X1hat = Xhat(1:M*N,indctr);
X2hat = Xhat(M*N+1:end,indctr);

axes(handles.axes1);
imagesc(Z);
colormap jet;
axis equal, axis off
title('Original Image')

smtDictSel = get(handles.smtDictSel,'SelectedObject');
smtDictSel = get(smtDictSel,'string');
D1_Type = strcmp(smtDictSel,'DFT');

axes(handles.axes2);
if D1_Type == 0
Tmp = Dmult(X1hat,M,N);
else
Tmp = Fmult(X1hat,M,N);    
end
imagesc( reshape(Tmp,M,N) ); colormap jet
axis equal, axis off
title('Recovered Smooth')

sprsDictSel = get(handles.sprsDictSel,'SelectedObject');
sprsDictSel = get(sprsDictSel,'string');
D2_Type = strcmp(sprsDictSel,'Marr Wavelet');

axes(handles.axes3);
if D2_Type == 0
    Tmp = X2hat;
else
    Tmp = Cmult(handles.dMtx, handles.maskMtx,X2hat,M, N);
end
imagesc( reshape(Tmp,M,N) ); 
     
colormap jet
axis equal, axis off
title('Recovered Sparse');



function displayBox_Callback(hObject, eventdata, handles)
% hObject    handle to displayBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of displayBox as text
%        str2double(get(hObject,'String')) returns contents of displayBox as a double


% --- Executes during object creation, after setting all properties.
function displayBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displayBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SptRsl_Callback(hObject, eventdata, handles)
% hObject    handle to SptRsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SptRsl as text
%        str2double(get(hObject,'String')) returns contents of SptRsl as a double


% --- Executes during object creation, after setting all properties.
function SptRsl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SptRsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when uipanel6 is resized.
function uipanel6_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
