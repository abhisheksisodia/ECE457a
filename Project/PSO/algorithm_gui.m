function varargout = algorithm_gui(varargin)
% ALGORITHM_GUI MATLAB code for algorithm_gui.fig
%      ALGORITHM_GUI, by itself, creates a new ALGORITHM_GUI or raises the existing
%      singleton*.
%
%      H = ALGORITHM_GUI returns the handle to a new ALGORITHM_GUI or the handle to
%      the existing singleton*.
%
%      ALGORITHM_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALGORITHM_GUI.M with the given input arguments.
%
%      ALGORITHM_GUI('Property','Value',...) creates a new ALGORITHM_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before algorithm_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to algorithm_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help algorithm_gui

% Last Modified by GUIDE v2.5 27-Jul-2015 23:02:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @algorithm_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @algorithm_gui_OutputFcn, ...
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


% --- Executes just before algorithm_gui is made visible.
function algorithm_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to algorithm_gui (see VARARGIN)

% Choose default command line output for algorithm_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes algorithm_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = algorithm_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function Lx_Callback(hObject, eventdata, handles)
% hObject    handle to Lx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lx as text
%        str2double(get(hObject,'String')) returns contents of Lx as a double


% --- Executes during object creation, after setting all properties.
function Lx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ly_Callback(hObject, eventdata, handles)
% hObject    handle to Ly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ly as text
%        str2double(get(hObject,'String')) returns contents of Ly as a double


% --- Executes during object creation, after setting all properties.
function Ly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function D_Callback(hObject, eventdata, handles)
% hObject    handle to D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D as text
%        str2double(get(hObject,'String')) returns contents of D as a double


% --- Executes during object creation, after setting all properties.
function D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pwt_Callback(hObject, eventdata, handles)
% hObject    handle to Pwt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pwt as text
%        str2double(get(hObject,'String')) returns contents of Pwt as a double


% --- Executes during object creation, after setting all properties.
function Pwt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pwt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the parameters from the GUI
Lx = str2num(get(handles.Lx,'String'));
Ly = str2num(get(handles.Ly,'String'));
D = str2num(get(handles.D,'String'));
Pwt = str2num(get(handles.Pwt,'String'));

%run the algorithm 10 times for the gui
xkrow = zeros(10,1);
xkcol = zeros(10,1);
N = zeros(10,1);
for i=1:10
    [xkrow(i,1),xkcol(i,1),N(i,1)] = pso(Lx,Ly,Pwt,D);
end
%output the results to the GUI
set(handles.krow,'String',mean(xkrow(:)));
set(handles.kcol,'String',mean(xkcol(:)));
set(handles.N,'String',mean(N(:)));

function krow_Callback(hObject, eventdata, handles)
% hObject    handle to krow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of krow as text
%        str2double(get(hObject,'String')) returns contents of krow as a double


% --- Executes during object creation, after setting all properties.
function krow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to krow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kcol_Callback(hObject, eventdata, handles)
% hObject    handle to kcol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kcol as text
%        str2double(get(hObject,'String')) returns contents of kcol as a double


% --- Executes during object creation, after setting all properties.
function kcol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kcol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N_Callback(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N as text
%        str2double(get(hObject,'String')) returns contents of N as a double


% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Determine the selected data set.
% str = get(hObject, 'String');
% val = get(hObject,'Value');
% % Set current data to the selected data set.
% switch str{val};
% case 'PSO' % User selects PSO.
%    set(handles.popupmenu1,'String',1);
% end
% % Save the handles structure.
% guidata(hObject,handles)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
