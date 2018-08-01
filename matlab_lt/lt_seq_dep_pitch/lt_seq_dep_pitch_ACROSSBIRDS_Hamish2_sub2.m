function [OUTSTRUCT] = ...
    lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub2(SeqDepPitch_AcrossBirds, ...
    birdtoget, expttoget, sylthis, PBSorMUSC, daylist, PCtimewindows, PCtimeWindowUsingWN)
%%
% PCtimewindows are hand designated windows
% NOTE, ffmat output is new ffmat recalcualted using hand designated
% windows (if they are defined);

%%

%% lt 7/14/18 - extracts pitch contours for desired days

birdname = SeqDepPitch_AcrossBirds.birds{birdtoget}.birdname;
exptname = SeqDepPitch_AcrossBirds.birds{birdtoget}.experiment{expttoget}.ExptID;

fdirthis = ['/bluejay5/lucas/analyses/seqdeppitch/GetPitchCont/' birdname '_' exptname];
% PBSorMUSC = 'PBS';


%%
i = birdtoget;
ii = expttoget;
DatThis = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning;
numdays = length(daylist);

% --------- TO COLLECT ACROSS DAYS
All_PCmat = cell(numdays,1);
All_twind = nan(numdays,2);
All_tbins = cell(numdays,1);
All_dayind = nan(numdays,1);
All_ffvals = cell(numdays,1);


for nn=1:numdays
    daythis = daylist(nn);
    
    
    % ----------- 1) load dat (PCmat
    fnamethis = [fdirthis '/syl' sylthis '_' PBSorMUSC '_day' num2str(daythis) '.mat'];
    PCmat = load(fnamethis);
    
    % ------------ 2) figure out time window for PC
    prms = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params;
    
    indtmp_sylthis = strcmp(DatThis.Params.PlotLearning.SylFieldsAll, sylthis);
    if length(DatThis.Params.PlotLearning.SylFieldsAll) == length(prms.SeqFilter.pc_time_window_list{daythis});
        % =================== version 1 -
        
        if isfield(prms, 'RecalculateFF')
            twind = prms.RecalculateFF.pc_time_window_list(:, indtmp_sylthis);
        else
            twind = prms.SeqFilter.pc_time_window_list{daythis}(:, indtmp_sylthis);
        end
        sylorigthis = prms.SeqFilter.OrigNoteID{daythis}(indtmp_sylthis);
        tbins = DatThis.Params.DayRawDat.pc_T{sylorigthis};
    else
        % ===== version 2 - I AM CONFUSED. GO BACK TO ORIGIANL SYL ID
        % i.e. find what origianl syl this dsired syl is, and use
        % originajmn t windows. these might not be the correct ones sicne
        % in some cases I recalcaulted FF using new windows.
        sylorig = prms.DayRawDat.syllables;
        if length(sylthis)==1
            syllowerthis = sylthis;
        else
            syllowerthis = lower(sylthis(regexp(sylthis, '[A-Z]'))); assert(length(syllowerthis)==1);
        end
        
        indorig = find(strcmp(sylorig, syllowerthis));
        twind = prms.DayRawDat.pc_time_window{indorig};
        if isfield(prms, 'RecalculateFF')
            disp('NOTE: recalc FF, the extracted twindow might not be accurate...');
        end
        tbins = DatThis.Params.DayRawDat.pc_T{indorig};
        
    end
    
    % ========= USE HAND LABELED WINDOWS?
    if exist('PCtimewindows', 'var')
        indtmp = [];
        
        if ~isempty(PCtimewindows)
            
            % ------------ FIRST, try using this syl in full (e.g. abB)
            indtmp_bird = find(strcmp(PCtimewindows, birdname));
            indtmp_expt = find(strcmp(PCtimewindows, exptname));
            indtmp_syl = find(strcmp(PCtimewindows, sylthis));
            
            indtmp = intersect(intersect(indtmp_bird+2, indtmp_syl), intersect(indtmp_expt+1, indtmp_syl));
            
            if PCtimeWindowUsingWN==1
                if isempty(indtmp)
                % then if indtmp is empty, it is because I did not enter
                % value for this syllablel. set time winodw to empty
                    twind = [];
                end
            else
                % -- then is empty must be becuase is coded as single syl
                % ...
                if isempty(indtmp)
                    % ---------- SECOND, if can't find, then use the lower syl
                    % (e.g. abB --> b)
                    assert(length(sylthis)>1, 'why? expect this to not be single syl, since failed to find match just now');
                    syllower = lower(sylthis(regexp(sylthis, '[A-Z]')));
                    indtmp_syl = find(strcmp(PCtimewindows, syllower));
                    
                    indtmp = intersect(intersect(indtmp_bird+2, indtmp_syl), intersect(indtmp_expt+1, indtmp_syl));
                end
                
                %             indtmp = intersect(find(strcmp(PCtimewindows, birdname))+1, find(strcmp(PCtimewindows, exptname)));
                
                assert(length(indtmp)==1, 'this bird, expt, syl combination doesnt exist');
                twindsec = PCtimewindows{indtmp+1};
                
                %         disp(['--- ' num2str(twindsec)]);
                % - convert from sec to samp
                if isempty(twindsec)
                    twind = [];
                else
                    tmp = find(tbins>=twindsec(1) & tbins<=twindsec(2));
                    twind = [tmp(1) tmp(end)];
                end
                disp('USING HAND CODED TWIND!!!');
            end
        end
    end
    
    
    % ################################# SAVE OUTPUT
    if strcmp(PBSorMUSC, 'MUSC')
        %             All_PCmat{nn} = PCmat.PCmat_tosave_MUSC(:, twind(1):twind(end));
        All_PCmat{nn} = PCmat.PCmat_tosave_MUSC;
    elseif strcmp(PBSorMUSC, 'PBS')
        %             All_PCmat{nn} = PCmat.PCmat_tosave(:, twind(1):twind(end));
        All_PCmat{nn} = PCmat.PCmat_tosave;
    end
    
        All_tbins{nn} = tbins;
        All_dayind(nn) = daythis;
    
    % ################################# SAVE OUTPUT THAT DEPENDS ON HAVING DEFINED TIME WINDOW
    if isempty(twind)
       
        % ---------- then output nans
        % don't do anything since already nans and empty.
%         All_twind(nn,:) = [];
%         All_ffvals{nn} = {};
        
    else
        % ============ confirm that FF correct
        if strcmp(PBSorMUSC, 'MUSC')
            fftmp = DatThis.AllDays_PlotLearning.DataMatrix_MUSC.(sylthis).FFvals_WithinTimeWindow{daythis};
            ffnew = mean(PCmat.PCmat_tosave_MUSC(:, twind(1):twind(end)),2);
        elseif strcmp(PBSorMUSC, 'PBS')
            fftmp = DatThis.AllDays_PlotLearning.DataMatrix.(sylthis).FFvals_WithinTimeWindow{daythis};
            ffnew = mean(PCmat.PCmat_tosave(:, twind(1):twind(end)),2);
        end
                
        All_twind(nn,:) = twind';
        All_ffvals{nn} = ffnew;
    end
    
    
    
    % ============================= DEBUGGING STUFF
    % -------- as an aside, FF for this data
    if (0)
        assert(strcmp(PBSorMUSC, 'PBS'));
        fftmp = DatThis.AllDays_PlotLearning.DataMatrix.(sylthis).FFvals_WithinTimeWindow{daythis};
        lt_figure; hold on;
        lt_plot_histogram(fftmp);
        assert(length(fftmp) == size(PCmat_tosave, 1));
    end
    
    % ---
    if (0)
        lt_figure; hold on;
        
        plot(tbins, PCmat.PCmat_tosave', '-k');
        tt = tbins(twind);
        line([tt(1) tt(1)], ylim, 'Color', 'r');
        line([tt(end) tt(end)], ylim, 'Color', 'r');
    end
    
end

%% == OUTSTRUCT
OUTSTRUCT.All_PCmat = All_PCmat;
OUTSTRUCT.All_twind = All_twind;
OUTSTRUCT.All_tbins = All_tbins;
OUTSTRUCT.All_dayind = All_dayind;
OUTSTRUCT.All_ffvals = All_ffvals;

