function [vals,fv]=evtaf_amp(bt,fbins,NT,NFFT,CS,USEX,USEX_TRIG);
% vals=evtaf_amp(bt,fbins,NT,NFFT,CS,USEX,USEX_TRIG);
% vals=[vals;fn2datenum(fn),tempvals,IND,filenum,TRIG];
%
%  bt - batch file
%  fbins - [Min Freq, Max Freq] for AMP comp
%  NT - target note
%  NFFT - length of template 
%  CS - chan spec
%  USEX - if == 1 look in X.rec for AMP calc times
% USEX_TRIG ==1 to use X.rec for which notes escaped ans which hit
% returns vals - 
%   vals=[datenum of file , FREQ vals , Note Index , file number,TRIG];
% usage example:
%   vals=evtaf_amp('batch.train',[5e2,1e4],'a',128,'obs1',1,0);
%   vals=evtaf_amp(bt,fbins,NT,NFFT,CS,TEMP1,TEMP2);
%   bt='batch.train';NT='a';TEMP1=1;TEMP2=0;NFFT=128;CS='obs1';fbins=[500,1e4];
%

fv=[];

fv=[];

NFFT=NFFT*2;
vals=[];
ff=load_batchf(bt);
for ii=1:length(ff)
	fn=ff(ii).name;
	ppp=findstr(fn,'.cbin');
	pppp=findstr(fn,'.');
	tmp=find(pppp<ppp);pppp=pppp(tmp(end));
	filenum=str2num(fn(pppp+1:ppp-1));
	if (~exist(fn,'file'))
		continue;
	end
	if (~exist([fn,'.not.mat'],'file'))
		continue;
	end
	rd=readrecf(fn,USEX);
	load([fn,'.not.mat']);

	if (USEX==USEX_TRIG)
		real_tt=rd.ttimes;
	else
		rd2=readrecf(fn,USEX_TRIG);
		real_tt=rd2.ttimes;
		clear rd2;
	end

	[dat,fs]=evsoundin('',fn,CS);

	if (~isfield(rd,'ttimes'))
		tt=[];
	else
		tt=rd.ttimes*1e-3;
	end

	for jj=1:length(tt)
		pp=find((onsets<=tt(jj)*1e3)&(offsets>=tt(jj)*1e3));
		if (~strcmp(labels(pp),NT))
			disp(['hey:',num2str(length(pp)),' ',labels(pp)]);

			continue;
		end
		IND=pp(1);

		ISTRIG=0;
		if (length(find((real_tt>=onsets(IND))&(real_tt<=offsets(IND))))>0)
			ISTRIG=1;
		end

		inds=fix(tt(jj)*fs)+[-(NFFT-1):0];
		dat2=dat(inds);

		fdat=(abs(fft(hamming(NFFT).*dat2))./NFFT);

		fv(length(fv)+1).fdat=fdat;

		ffv=get_fft_freqs(NFFT,fs);
		ffv=ffv(1:end/2);
		inds2=find((ffv>=fbins(1))&(ffv<=fbins(2)));

		vals=[vals;fn2datenum(fn),sum(fdat(inds2)),IND,filenum,ISTRIG];
	end
end
return;
