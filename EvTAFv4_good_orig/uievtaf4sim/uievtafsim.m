function varargout = uievtafsim(varargin)
% 
%   RUN AS uievtafsim(batchfile,configfile);

%      UIEVTAFSIM, by itself, creates a new UIEVTAFSIM or raises the existing
%      singleton*.
%
%      H = UIEVTAFSIM returns the handle to a new UIEVTAFSIM or the handle to
%      the existing singleton*.
%
%      UIEVTAFSIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UIEVTAFSIM.M with the given input arguments.
%
%      UIEVTAFSIM('Property','Value',...) creates a new UIEVTAFSIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before uievtafsim_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to uievtafsim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help uievtafsim

% Last Modified by GUIDE v2.5 18-Aug-2009 11:42:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @uievtafsim_OpeningFcn, ...
    'gui_OutputFcn',  @uievtafsim_OutputFcn, ...
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
return;

% --- Executes just before uievtafsim is made visible.
function uievtafsim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to uievtafsim (see VARARGIN)

% Choose default command line output for uievtafsim
handles.output = hObject;

% load up the files in the batch file - ignore files that do not exist
fid=fopen(varargin{1},'r');
if (fid==-1)
    disp(['Bad batchfile']);
    return;
end

ind=0;
while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    ind=ind+1;
    files{ind}=fn;
end
if (ind <1)
    disp(['No input files found']);
    return;
end

NDFile=varargin{2};
[ND,OP]=ReadEvTAFv4ConfigFile(NDFile,0);

for ijk=1:length(ND)
    if (size(ND(ijk).Templ,1)==0)
        disp(['If you have the templates in a file run '])
        disp(['AddTemplatesToEvConfig(''config.evconfig2'',''Template.dat''']);
        disp([' or ']);
        disp(['Size of the templates in this config file is 0']);
        disp(['You need to run EvTAF (meaning hit the GO button)']);
        disp(['Then stop EvTAF (hit the STOP button']);
        disp(['Restart EvTAF (the run arrow on the top LabVIEW window bar']);
        disp(['and Save the config file again to get the templates written into the config file']);
    end
end

NDIndex=1;
handles.NDIndex=1;
handles.ND=ND;
handles.OP=OP;
handles.t=0;handles.fs=0;handles.dat=0;handles.sp=[];handles.labels='';
handles.onsets=[];handles.offsets=[];
guidata(hObject, handles);


tmpstr=[];
for ijk=1:length(handles.ND)
    tmpstr{ijk}=num2str(ijk);
end
set(handles.NDIndexBox,'String',tmpstr);
set(handles.NDIndexBox,'Value',handles.NDIndex);
guidata(hObject, handles);

[NDOut,NDIndexOut,WasCanceled]=CounterSetup(handles.ND,handles.NDIndex);
if WasCanceled
    return;
end
handles.ND=NDOut;
handles.NDIndex=NDIndexOut;

tmpstr=[];
for ijk=1:length(handles.ND)
    tmpstr{ijk}=num2str(ijk);
end
set(handles.NDIndexBox,'String',tmpstr);
set(handles.NDIndexBox,'Value',handles.NDIndex);
guidata(hObject, handles);

templates = ND(NDIndex).Templ;
handles.NTempl = size(templates,2);
handles.TemplLen = size(templates,1);

handles.files=files;
handles.Nfile=1;
handles.CHANSPEC=get(handles.ChanSpecBox,'String');
handles.SEGTH=str2num(get(handles.SegThreshBox,'String'));
set(handles.SPMinLevel,'Value',0.7);
set(handles.SPMaxLevel,'Value',1);
set(handles.RefracPerBox,'String',num2str(handles.ND(NDIndex).TrigRefrac));

handles.TrigT=[];
handles.TRIGPLT=-1;
guidata(hObject, handles);

UIEvTAFv4Sim_PlotData(hObject,handles,1);
handles=guidata(hObject);

axes(handles.ValAxes);ylim([0,5]);
linkaxes([handles.SpecGramAxes,handles.LabelAxes,handles.ValAxes],'x');
% UIWAIT makes uievtafsim wait for user response (see UIRESUME)
uiwait(handles.uievtafsim);
return;

% --- Outputs from this function are returned to the command line.
function varargout = uievtafsim_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = [];
delete(hObject);
return;

% --- Executes on button press in PrevFileBtn.
function PrevFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to PrevFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.Nfile>1)
    if (get(handles.WriteRecFileToggle,'Value')==get(handles.WriteRecFileToggle,'Max'))
        DoWriteRecFile(hObject, handles);
    end
    handles.Nfile=handles.Nfile-1;
    UIEvTAFv4Sim_PlotData(hObject,handles,1);
end
return;

% --- Executes on button press in NextFileBtn.
function NextFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to NextFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.Nfile<length(handles.files))
    if (get(handles.WriteRecFileToggle,'Value')==get(handles.WriteRecFileToggle,'Max'))
        DoWriteRecFile(hObject, handles);
    end
    handles.Nfile=handles.Nfile+1;
    UIEvTAFv4Sim_PlotData(hObject,handles,1);
end
return;

% --- Executes on button press in MoveRightBtn.
function MoveRightBtn_Callback(hObject, eventdata, handles)
% hObject    handle to MoveRightBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SpecGramAxes);
v=axis;
xlim(v(1:2)-0.9*(v(2)-v(1)));
return;

% --- Executes on button press in MoveLeftBtn.
function MoveLeftBtn_Callback(hObject, eventdata, handles)
% hObject    handle to MoveLeftBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SpecGramAxes);
v=axis;
xlim(v(1:2)+0.9*(v(2)-v(1)));
return;

% --- Executes on button press in SkipToFileBtn.
function SkipToFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SkipToFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[nf,canc]=SkipToFile(handles.files,handles.Nfile);
if (canc==0)
    if (get(handles.WriteRecFileToggle,'Value')==get(handles.WriteRecFileToggle,'Max'))
        DoWriteRecFile(hObject, handles);
    end
    handles.Nfile=nf;
    UIEvTAFv4Sim_PlotData(hObject,handles,1);
end
return;

function SegThreshBox_Callback(hObject, eventdata, handles)
% hObject    handle to SegThreshBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SegThreshBox as text
%        str2double(get(hObject,'String')) returns contents of SegThreshBox as a double

handles.SEGTH=str2num(get(handles.SegThreshBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function SegThreshBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SegThreshBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

% --- Executes on button press in ZoomOnBtn.
function ZoomOnBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomOnBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom on;
return;

% --- Executes on button press in PanOnBtn.
function PanOnBtn_Callback(hObject, eventdata, handles)
% hObject    handle to PanOnBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pan xon;
return;

function ChanSpecBox_Callback(hObject, eventdata, handles)
% hObject    handle to ChanSpecBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ChanSpecBox as text
%        str2double(get(hObject,'String')) returns contents of ChanSpecBox
%        as a double

handles.CHANSPEC=get(handles.ChanSpecBox,'String');
guidata(hObject,handles);
UIEvTAFv4Sim_PlotData(hObject,handles,1);
return;

% --- Executes during object creation, after setting all properties.
function ChanSpecBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChanSpecBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;


% --- Executes on button press in ZoomXonBtn.
function ZoomXonBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomXonBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom off;
zoom xon;
return;

% --- Executes on button press in ReplotBtn.
function ReplotBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ReplotBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UIEvTAFv4Sim_PlotData(hObject,handles,1);
return;

% --- Executes during object creation, after setting all properties.
function FileNumBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileNumBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;


% --- Executes on button press in QuitBtn.
function QuitBtn_Callback(hObject, eventdata, handles)
% hObject    handle to QuitBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.WriteRecFileToggle,'Value')==get(handles.WriteRecFileToggle,'Max'))
    DoWriteRecFile(hObject, handles);
end
uiresume(handles.uievtafsim);
return;

% --- Executes on button press in UseSim.
function UseSim_Callback(hObject, eventdata, handles)
% hObject    handle to UseSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseSim

PlotTafVals(hObject,handles);

return;


% --- Executes on button press in SetCntrBtn.
function SetCntrBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SetCntrBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[NDOut,NDIndexOut,WASCANC]=CounterSetup(handles.ND,handles.NDIndex);

if (~WASCANC)
    handles.ND=NDOut;
    handles.NDIndex=NDIndexOut;
    guidata(hObject, handles);
    tmpstr=[];
    for ijk=1:length(handles.ND)
        tmpstr{ijk}=num2str(ijk);
    end
    set(handles.NDIndexBox,'String',tmpstr);
    set(handles.NDIndexBox,'Value',handles.NDIndex);
    guidata(hObject, handles);
    
    UIEvTAFv4Sim_PlotData(hObject,handles,1);
end
return;

% --- Executes on button press in UseXTmpFile.
function UseXTmpFile_Callback(hObject, eventdata, handles)
% hObject    handle to UseXTmpFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseXTmpFile

if (get(hObject,'Value')==get(hObject,'Max'))
    handles.USEXTMP=1;
    guidata(hObject, handles);
else
    handles.USEXTMP=0;
    guidata(hObject, handles);
end
PlotTafVals(hObject,handles);
return;

% --- Executes on slider movement.
function SPMinLevel_Callback(hObject, eventdata, handles)
% hObject    handle to SPMinLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if (get(handles.SPMinLevel,'Value')<get(handles.SPMaxLevel,'Value'))
    axes(handles.SpecGramAxes);
    temp1=get(handles.SPMinLevel,'Value');temp2=get(handles.SPMaxLevel,'Value');
    mn=handles.SPMin;mx=handles.SPMax;
    caxis([temp1*(mx-mn)+mn,temp2*mx]);
else
    set(handles.SPMinLevel,'Value',.999*get(handles.SPMaxLevel,'Value'));
end
return;

% --- Executes during object creation, after setting all properties.
function SPMinLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SPMinLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
return;

% --- Executes on slider movement.
function SPMaxLevel_Callback(hObject, eventdata, handles)
% hObject    handle to SPMaxLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if (get(handles.SPMaxLevel,'Value')>get(handles.SPMinLevel,'Value'))
    axes(handles.SpecGramAxes);
    temp1=get(handles.SPMinLevel,'Value');temp2=get(handles.SPMaxLevel,'Value');
    mn=handles.SPMin;mx=handles.SPMax;
    caxis([temp1*(mx-mn)+mn,temp2*mx]);
else
    set(handles.SPMaxLevel,'Value',1.001*get(handles.SPMinLevel,'Value'));
end
return;

% --- Executes during object creation, after setting all properties.
function SPMaxLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SPMaxLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
return;

function RefracPerBox_Callback(hObject, eventdata, handles)
% hObject    handle to RefracPerBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RefracPerBox as text
%        str2double(get(hObject,'String')) returns contents of RefracPerBox as a double

handles.REFRAC=str2num(get(handles.RefracPerBox,'String'));
guidata(hObject, handles);
PlotTafVals(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function RefracPerBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RefracPerBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

% --- Executes on selection change in NDIndexBox.
function NDIndexBox_Callback(hObject, eventdata, handles)
% hObject    handle to NDIndexBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns NDIndexBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NDIndexBox

NDIndex=get(handles.NDIndexBox,'Value');
handles.NDIndex=NDIndex;
guidata(hObject, handles);

%NDOut=CounterSetup(handles.ND,handles.NDIndex);
%handles.ND=NDOut;
%guidata(hObject, handles);

ND=handles.ND;
NDIndex=handles.NDIndex;

templates = ND(NDIndex).Templ;
handles.NTempl = size(templates,2);
handles.TemplLen = size(templates,1);
set(handles.RefracPerBox,'String',num2str(handles.ND(NDIndex).TrigRefrac));
handles.TrigT=[];
handles.TRIGPLT=-1;
set(handles.RefracPerBox,'String',num2str(handles.ND(handles.NDIndex).TrigRefrac));
guidata(hObject, handles);

UIEvTAFv4Sim_PlotData(hObject,handles,0);

return;

% --- Executes during object creation, after setting all properties.
function NDIndexBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NDIndexBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TriggerListBox.
function TriggerListBox_Callback(hObject, eventdata, handles)
% hObject    handle to TriggerListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns TriggerListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TriggerListBox


% --- Executes during object creation, after setting all properties.
function TriggerListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TriggerListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in WriteEvConfigButton.
function WriteEvConfigButton_Callback(hObject, eventdata, handles)
% hObject    handle to WriteEvConfigButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[configfn,pth]=uiputfile('*.evconfig2','Save evconfig2 file');
configfn=fullfile(pth,configfn);
WriteEvTAFv4ConfigFile(configfn,handles.ND,handles.OP);
disp(['Config file written to : ',configfn]);
return

function DoWriteRecFile(hObject, handles)
fn=handles.files{handles.Nfile};
ND=handles.ND;
rd=readrecf(fn,0);
SimTrigInfo=[];
for iND=1:length(ND)
    SimTrigInfo=[SimTrigInfo;...
        ND(iND).TriggerTimes*1e3,ones(size(ND(iND).TriggerTimes))*(iND-1)];
end
[sortv,sorti]=sort(SimTrigInfo(:,1));
SimTrigInfo=SimTrigInfo(sorti,:);

rd.ttimes=SimTrigInfo(:,1);
rd.trignote=SimTrigInfo(:,2);

outf=wrtrecf(fn,rd,1);
disp(['Rec file written to :',outf]);
return;


% --- Executes on button press in WriteRecFileToggle.
function WriteRecFileToggle_Callback(hObject, eventdata, handles)
% hObject    handle to WriteRecFileToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of WriteRecFileToggle


