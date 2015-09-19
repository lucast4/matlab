function ND=EvTAFv4Sim(batchf,configfile,ChanSpec);
%
%

[ND,OP]=ReadEvTAFv4ConfigFile(configfile);

% do some checks to make sure everything is ok
for ii=1:length(ND)
    if isempty(ND(ii).Templ)
        disp(['Templates are empty for note detect info index = ',num2str(ii)]);
        if exist(ND(ii).TemplFile,'file')
            disp(['Load the templates from the template file : ',ND(ii).TemplFile]);
            ND(ii).Templ=load(ND(ii).TemplFile);
        else
            disp(['Could not file the Template File indicated : ',ND(ii).TemplFile]);
            disp(['Probably a path mismatch']);
            disp(['Please choose the right template file']);
            [FileName,PathName,FilterIndex] = uigetfile('*',['Template for Note Index ',num2str(ii)]);
            ND(ii).TemplFile=[FileName,PathName];
            ND(ii).Templ=load(ND(ii).TemplFile);
        end
    end

    %verify template normalizations
    for jj=1:size(ND(ii).Templ,2)
        ND(ii).Templ(1:6,jj)=0;
        ND(ii).Templ(:,jj)=ND(ii).Templ(:,jj)./max(ND(ii).Templ(:,jj));
    end

    % get the number of data points that have to be read in based on the
    % size of the templates
    if (ii==1)
        NFFT=2*size(ND(ii).Templ,1);
    else
        if (NFFT~=2*size(ND(ii).Templ,1))
            disp(['ERROR - EvTAF Limitation all templates have the be the same length in freq bins']);
            return;
        end
    end

    if (length(ND(ii).CntRng)~=size(ND(ii).Templ,2))
        disp(['Counter Range values and Templates are of different sizes for Note index : ',num2str(ii)]);
        disp(['Cannot continue']);
        return;
    end
end


inputfiles=[];
fid=fopen(batchf,'r');
while (1)
    fn=fgetl(fid);
    if (ischar(fn))
        if (exist(fn,'file'))
            inputfiles(length(inputfiles)+1).fn=fn;
        else
            disp(['Could not fine file : ',fn]);
            disp('Skipping it');
        end
    else
        break;
    end
end
fclose(fid);

for IFile=1:length(inputfiles)
    %loop over every file in the batchfile
    fn=inputfiles(IFile).fn;
    disp(fn);
    [dat,fs]=evsoundin('',fn,ChanSpec); % load up the raw data

    PREBUFIND=ceil(OP.FileBufferLeng*fs);

    %skip the data in the prebuffer - Cannot trigger off this data since
    %it is just held in the buffer and no sound playback happens if the box
    %is 'quiet'

    % initialize the counters
    for iND=1:length(ND)
        for iCnt=1:length(ND(iND).CntRng)
            if (ND(iND).CntRng(iCnt).Mode==1)
                %this is evtaf mode it means that the counters count up as
                %long as the template comparison is below the threshold
                %CntRng(iCnt).TH
                % this mode keeps track of consequtive template matches

                %init counters such they are in a state which should not
                %immeidately trigger
                ND(iND).Cnt(iCnt)=0;
            else
                % BirdTAF mode, counters go to zero and stay at zero as
                % long as template comparison is below threhold and count
                % up as long as it is above match treshold
                % this mode keeps track of how long ago a template match
                % occurred
                % one addition - BTMin determines minimum Counter value
                % requred before the counters can clear again
                % set to zero to have normal birdtaf function set higher to
                % avoid some spurious triggers caused by noise

                ND(iND).Cnt(iCnt)=ND(iND).CntRng(iCnt).Max+1;
            end
            ND(iND).Detect=0; % was the note detected
            ND(iND).Trigger=0; % note detected, amp and ff ranges match -> get trigger
            ND(iND).TimeSinceDet=ND(iND).TrigRefrac+1; % can immediately trigger since they clear refrac requirement
            ND(iND).TimeSinceTrig=ND(iND).TrigRefrac+1;
            ND(iND).TriggerTimes=[];
        end
    end

    DeltaT=NFFT/fs; % time length of each data frame in seconds
    %load up the data
    for dataind=PREBUFIND:NFFT:(length(dat)-NFFT+1)
        CurrentTime=(dataind+NFFT-1)/fs;
        % take data in chunks 2x the length of the templates
        tempdat=dat(dataind:dataind+NFFT-1);

        % FFT the data chunk
        fftdat=abs(fft(tempdat.*hamming(length(tempdat))));
        fftdat=fftdat(1:fix(NFFT/2));

        %normalize to remove low freq noise and make max point = 1.0
        fftdat(1:6)=0;fftdat=fftdat./max(fftdat);

        % go through each template set in Note Detection Info
        %(see EvTAFv4 LabVIEW window
        for iND=1:length(ND)
            NTemplates=size(ND(iND).Templ,2);

            MatchArray=zeros([1,NTemplates]);
            for itempl=1:NTemplates
                Template=ND(iND).Templ(:,itempl);
                MINV=ND(iND).CntRng(itempl).Min;
                MAXV=ND(iND).CntRng(itempl).Max;
                MODE=ND(iND).CntRng(itempl).Mode;
                DONOT=ND(iND).CntRng(itempl).Not;
                BTMIN=ND(iND).CntRng(itempl).BTAFmin;
                TH=ND(iND).CntRng(itempl).TH;
                CntVal=ND(iND).Cnt(itempl);

                % template comparison
                comparval = sum((fftdat - Template).^2);

                %update the counters
                if (comparval <= TH)
                    %template 'match'
                    if (MODE==1)
                        % EvTAF Mode means count up
                        CntVal=CntVal+1;
                    elseif (MODE==0)
                        %BirdTAF mode means zero out
                        if (CntVal>=BTMIN)
                            % when using long counts between templates this
                            % can help avoid some false resets of the bird
                            % taf trigger.  in essence this sauys the the
                            % trigger had to happen more than BTMIN counts
                            % ago in order to clear the counter now
                            % SET TO ZERO TO HAVE USUAL BIRDTAF BEHAVIOR
                            CntVal=0;
                        else
                            CntVal=CntVal+1;
                        end
                    end
                else
                    %template 'mismatch'
                    if (MODE==1)
                        % EvTAF Mode means zero out
                        CntVal=0;
                    elseif (MODE==0)
                        %BirdTAF mode means count up
                        CntVal=CntVal+1;
                    end
                end
                ND(iND).Cnt(itempl)=CntVal;
                %end update counters

                %check if this counter is in range - save the value in
                %MatchArray
                if ((CntVal>=MINV)&(CntVal<=MAXV))
                    MatchArray(itempl)=1;
                else
                    MatchArray(itempl)=0;
                end

                %if this is a NOT template then match happens when counters
                %are NOT is this range so invert the answer above
                if (DONOT==1)
                    MatchArray(itempl) = ~MatchArray(itempl);
                end
            end

            % now do the abitrary logic function listed in 'Counter
            % Logic' in the labview window
            % this may get a little strange - user is not supposed to
            % change the field 'VarName' in 'Counter Rng' but i put
            % this in to make sure that we're covered if they do
            % do this part in a seperate function to aviod any same
            % named variables in the program
            DidItMatchLogic=CheckLogicFunc(ND(iND),MatchArray);

            %If it matches then set Detect to 1 unless refrac violation or
            %note already detected (detect already = 1)
            if (DidItMatchLogic)
                if (ND(iND).TimeSinceDet>ND(iND).TrigRefrac)
                    if (~ND(iND).Detect)
                        ND(iND).Detect=1;
                        %if detect was already 1 then do not reset detect
                        %timer
                        ND(iND).TimeSinceDet=0;
                    end
                end
            end

            % if note was detected then check to see if it's Delay time
            % critereon is met yet?
            %
            if (ND(iND).Detect==1)
                if (ND(iND).TimeSinceDet >= ND(iND).DelayToContingen)
                    % Check to see if Fund Freq passes contingency
                    MinFreq=ND(iND).FreqRng.MinFreq;
                    MaxFreq=ND(iND).FreqRng.MaxFreq;
                    FreqTHMin=ND(iND).FreqRng.FreqTHMin;
                    FreqTHMax=ND(iND).FreqRng.FreqTHMax;
                    NBins=ND(iND).FreqRng.NBins;
                    FreqMode=ND(iND).FreqRng.FreqMode;

                    %frequency bins corresponding to the FFT data
                    FreqVals=get_freqvals(NFFT,fs);

                    %get the indicies for the FFT that are between MinFreq
                    %and MaxFreq for the Freq Bounds (i know this is a
                    %silly way to do this but i'm trying to replicate the
                    %labview so that this is as accurate a representation
                    %as possible
                    tempind=1;
                    while (1)
                        if (FreqVals(tempind)>=MinFreq)
                            break;
                        end
                        tempind=tempind+1;
                    end
                    MinInd=tempind;

                    tempind=1;
                    while (1)
                        if (FreqVals(tempind)>=MaxFreq)
                            break;
                        end
                        tempind=tempind+1;
                    end
                    MaxInd=tempind;
                    
                    % pull out the FFT dat in this frequency range
                    tempdat=fftdat(MinInd:MaxInd);
                    %find the maximum value
                    [maxv,maxi]=max(tempdat);
                    maxi=maxi+MinInd-1;

                    %take the Nbins on each side
                    useinds = maxi + [0:2*NBins] - NBins; % again matching labview

                    % i think labview does this for me, i did not do it
                    % explicitly there have to do it here in case
                    useinds=useinds(find((useinds>0)&(useinds<=length(fftdat))));

                    % get the normalized - smoothed peak frequency
                    FundFreq=sum(fftdat(useinds).*FreqVals(useinds).')./sum(fftdat(useinds));
                    %disp([num2str(iND),' FF = ',num2str(FundFreq)]);
                    %%%%%%%%%%%  AMP CONTINGECY %%%%%%%%%%%%%%%%%%
                    % for programming simplicity here i will do the volume
                    % (amp) detection here
                    MinAmpFreq=ND(iND).AmpRng.MinFreq;
                    MaxAmpFreq=ND(iND).AmpRng.MaxFreq;
                    AmpMode=ND(iND).AmpRng.AmpMode;
                    AmpTH=ND(iND).AmpRng.AmpThresh;
                    tempind=1;
                    while (1)
                        %used > NOT >= in LV here don't know why?
                        if (FreqVals(tempind)>MinAmpFreq)
                            break;
                        end
                        tempind=tempind+1;
                    end
                    MinInd=tempind;

                    tempind=1;
                    while (1)
                        if (FreqVals(tempind)>MaxAmpFreq)
                            break;
                        end
                        tempind=tempind+1;
                    end
                    MaxInd=tempind;

                    tempdat=fftdat(MinInd:MaxInd);
                    %sum up the power in the freq range to get vol.
                    Vol=sum(tempdat);
                    if (AmpMode==0)
                        %below hits
                        if (Vol <= AmpTH)
                            TRIGGER=1;
                        else
                            ND(iND).Detect=0;
                            TRIGGER=0;
                        end
                    else
                        %above hits
                        if (Vol >= AmpTH)
                            TRIGGER=1;
                        else
                            ND(iND).Detect=0;
                            TRIGGER=0;
                        end
                    end


                    %now check the FundFreq Range if Amp passed
                    if (FreqMode==0)
                        % below FreqTHMax hits
                        % if you ever get confused for this and the next
                        % mode make FreqTHMax and FreqThMin the same value
                        if (FundFreq<=FreqTHMax)
                            TRIGGER = TRIGGER & 1;
                        else
                            ND(iND).Detect=0;
                            TRIGGER=0;
                        end
                    elseif (FreqMode==1)
                        % above FreqTHMin hits
                        if (FundFreq>=FreqTHMin)
                            TRIGGER = TRIGGER & 1;
                        else
                            ND(iND).Detect=0;
                            TRIGGER=0;
                        end
                    elseif (FreqMode==2)
                        %inside range hits
                        if ((FundFreq>=FreqTHMin)&(FundFreq<=FreqTHMax))
                            TRIGGER = TRIGGER & 1;
                        else
                            ND(iND).Detect=0;
                            TRIGGER=0;
                        end
                    else (FreqMode==3)
                        % outside range hits
                        if (~((FundFreq>=FreqTHMin)&(FundFreq<=FreqTHMax)))
                            TRIGGER = TRIGGER & 1;
                        else
                            ND(iND).Detect=0;
                            TRIGGER=0;
                        end
                    end

                    if (TRIGGER==0)
                        ND(iND).Detect=0;
                    else
                        if (iND==1)
                            disp(['FF=',num2str(FundFreq)]);
			    disp([' cnt=',num2str(ND(iND).Cnt)]);
                        end
                        
                        if (ND(iND).TimeSinceTrig>ND(iND).TrigRefrac)
                            ND(iND).TimeSinceTrig=0;
                            ND(iND).TriggerTimes=[ND(iND).TriggerTimes;CurrentTime];
                            ND(iND).Detect=0;
                        end
                    end

                end
            end
        end % loop over template sets (ND)

        %incremement the timers for next loop iteration
        for iND=1:length(ND)
            ND(iND).TimeSinceDet=ND(iND).TimeSinceDet + DeltaT;
            ND(iND).TimeSinceTrig=ND(iND).TimeSinceTrig + DeltaT;
        end

    end %loop over data in this file

    for iND=1:length(ND)
        %disp(['IND  = ',num2str(iND)]);
        %for ijk=1:length(ND(iND).TriggerTimes)
        %    disp(num2str(ND(iND).TriggerTimes(ijk),6));
        %end
        disp([' ']);
        rd=readrecf(fn);
        tt=[rd.ttimes,zeros(size(rd.ttimes))];
        for ijk=1:size(tt,1)
            pppp=findstr(rd.pbname{ijk},'Templ = ');
            strtemp=rd.pbname{ijk};
            tt(ijk,2)=str2num(strtemp(pppp+8:end));
        end
        
        
        %compare the sim trigs to the actual trigs
	if (isempty(tt))
            disp(['No trigs in recfile : Nact = 0  - Nsim = ',num2str(length(ND(iND).TriggerTimes))]);
		continue;
	end
		
        pppp=find(tt(:,2)==iND-1);
        if (length(pppp)==0)
            continue;
        end
        tttemp=tt(pppp,1);
        if (length(tttemp)~=length(ND(iND).TriggerTimes))
            disp(['XTrig Times not the same size : Nact = ',num2str(length(tttemp)),' Nsim = ',num2str(length(ND(iND).TriggerTimes))]);
        end
        for ijk=1:length(ND(iND).TriggerTimes)
            [yy,iy]=min(abs(ND(iND).TriggerTimes(ijk)*1e3-tttemp));
            disp(['iND = ',num2str(iND)]);
            disp(['Act Trig Time = ',num2str(tttemp(iy)),' Sim Trigger Time = ',num2str(ND(iND).TriggerTimes(ijk)*1e3)]);
            disp(['Diff in ms = ',num2str( tttemp(iy) - ND(iND).TriggerTimes(ijk)*1e3)]);
        end
    end
end % loop over files




