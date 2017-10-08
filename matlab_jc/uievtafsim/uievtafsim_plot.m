function uievtafsim_plotdata(hObject,handles);
%uievtafsim_plotdata(hObject,handles);

set(handles.FileNumBox,'String',[num2str(handles.Nfile),'/',num2str(length(handles.files))]);

fn=handles.files{handles.Nfile};
[dat,fs]=evsoundin('',fn,handles.CHANSPEC);
[sm,sp,t,f]=evsmooth(dat,fs,0);
sp=abs(sp);

if (~exist([fn,'.not.mat'],'file'))
    min_int = 5;
    min_dur = 30;
    sm_win=2;
    [onsets,offsets]=SegmentNotesJC(sm,fs,min_int,min_dur,handles.SEGTH);
    onsets = onsets*1e3;offsets=offsets*1e3;
    labels = [];
else
    load([fn,'.not.mat']);
    %if (handles.SEGTH~=threshold)
    %    [onsets,offsets]=SegmentNotes(sm,Fs,min_int,min_dur,threshold);
    %    labels = char(ones([1,length(onsets)])*48);
    %end
end

axes(handles.SpecGramAxes);hold off;
imagesc(t,f,log(abs(sp)));
set(gca,'YD','n');
title(fn,'Interpreter','none');xlim([0,t(end)]);ylim([0,1e4]);

handles.SPMax=max(max(log(sp)));
handles.SPMin=min(min(log(sp)));
temp1=get(handles.SPMinLevel,'Value');temp2=get(handles.SPMaxLevel,'Value');
mn=handles.SPMin;mx=handles.SPMax;
caxis([temp1*(mx-mn)+mn,temp2*mx]);

handles.t=t;
handles.fs=fs;
handles.dat=dat;
guidata(hObject,handles);

PlotTafVals(hObject,handles);
handles=guidata(hObject);

axes(handles.LabelAxes);cla;hold off;
guidata(hObject,handles);
if (exist('handles.ActTrigTimes'))
    delete(handles.ActTrigTimes);
    handles.ActTrigTimes=[];
end
set(handles.LabelAxes,'XTick',[],'YTick',[]);

if (length(labels)>0)
    handles.LabelHandl=text((onsets+offsets).'*5e-4,0*onsets.',labels.');hold on;
else
    handles.LabelHandl=[];
end

rdata=readrecf(fn);
if (exist('rdata'))
    if ~isempty(rdata)
    if (length(rdata.ttimes)>0)
        handles.ActTrigTimes=plot(rdata.ttimes*1e-3,0*rdata.ttimes,'r^');hold on;
    end
    end
end

trigtimes = handles.TrigT;
tmp_plt=plot(trigtimes,0*trigtimes-0.5,'bs');
axis([0,t(end),-1,1]);
handles.TRIGPLT=tmp_plt;

linkaxes([handles.SpecGramAxes,handles.ValAxes,handles.LabelAxes],'x');

handles.tlim=[0,t(end)];
handles.sp=sp;
handles.fs=fs;
handles.onsets=onsets;
handles.offsets=offsets;
handles.TrigT=trigtimes;
guidata(hObject,handles);

return;
