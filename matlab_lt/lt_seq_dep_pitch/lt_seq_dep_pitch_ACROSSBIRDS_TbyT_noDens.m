function [NtargCell, NnontargCell, FFbinnedCell, TbinnedCell, FFsinglebinMat, ...
    Nratio_hilo_targ, NumDatPerRend, BirdNumMat] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_noDens(DATBYREND, ...
    dotrain, dosametype, onlyifSigLearn, densitymethod, xedges, mintime_fromref, ...
    maxtime_fromref, minrends_inbin, minsongs_inbin)
%% lt 7/8/18 - modified from lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S3
% HERE: doesnt' separate data by hi or lo density.


%%
maxnumsyls = max(DATBYREND.Sylcounter);
xcenters = xedges(1:end-1)+diff(xedges)/2;

% binsize = xedges(2)-xedges(1);
% xcenters = xedges(1:end-1)+binsize/2;

%%

% ================ OUTPUT
NtargCell = cell(maxnumsyls,1); % syls x density(0,1);
NnontargCell = cell(maxnumsyls,1);

FFbinnedCell = cell(maxnumsyls,1);
TbinnedCell = cell(maxnumsyls,1);
FFsinglebinMat = nan(maxnumsyls, 1);
NumDatPerRend = nan(maxnumsyls,1);
BirdNumMat = nan(maxnumsyls,1);

% lt_figure; hold on;
% ================ GO THRU EACH EXPT, COLLECT BINNED FF DEV AND DENSITY
for ss=1:maxnumsyls
    
        indsthis = DATBYREND.Sylcounter==ss & DATBYREND.IsDurTrain==dotrain ...
            & DATBYREND.IsSame==dosametype;
        
        if ~any(indsthis)
            continue
        end
        
        birdthis = unique(DATBYREND.Birdnum(indsthis)); assert(length(birdthis)==1);
        
        
        % ======================== COLLECT DENSITY OF TARG AND NONTARG
        if strcmp(densitymethod, 'default')
            n_targ = DATBYREND.Density_targ(indsthis);
            n_nontarg = DATBYREND.Density_nontarg(indsthis);
        elseif strcmp(densitymethod, 'FFwindow')
            disp('NOT FINISHED CODING YET ...');
            pause;
            
            % --- collect the number of rends of targ and nontarg in following window
            tdev = DATBYREND.Time_dev(indsthis);
            tdev_targ = DATBYREND.Time_dev_targ(indsthis);
            nrends = length(tdev);
            for r =1:nrends
                
                t = tdev{r}*(24*60);
                indstmp =t>-0.25 & t<0.25;
            end
        end
        

        
        % ======================== COLLECT BINNED FF DEV ACROSS ALL TRIALS
        tdev = cell2mat(DATBYREND.Time_dev(indsthis));
        tdev = tdev*(24*60); % convert from day to minutes
        ffdev = cell2mat(DATBYREND.FF_dev(indsthis));
        
        xbins = discretize(tdev, xedges);
        [ffbinned_mean, ffbinned_sem] = grpstats(ffdev, xbins, {'mean', 'sem'});
        tbinned = xcenters(unique(xbins));
        
        
        % =========================== POOL ALL DATA WITHIN THE DESIRED BIN
        if (0) % indirect, takes data from already made bins
            binsthis = find(xcenters>mintime_fromref & xcenters<maxtime_fromref);
            indsinbin = ismember(xbins, binsthis); % renditions that are in desired windwo
            
            ffdev_singlebin = mean(ffdev(indsinbin));
        else
            indstmp = tdev>mintime_fromref & tdev<maxtime_fromref;
            ffdev_singlebin = mean(ffdev(indstmp));
            
        end
        
        
        
        
        % ======================== CRITERIA TO PASS TO COLLECT DATA
        % --------- 1) Min rends in high and low density bins
        if (1) % method 1: bin number of cases in the bins within timeframe
            binsthis = find(xcenters>mintime_fromref & xcenters<maxtime_fromref);
            indsinbin = ismember(xbins, binsthis); % renditions that are in desired windwo
            
            % -------- only continue if enough renditions in bin
            nrendsinbin = sum(indsinbin);
            if nrendsinbin<minrends_inbin
                disp('skip - not enough rends in bin')
                continue
            end
        end
        
        
        % ----------- COLLECT AVERAGE NUMBER OF DATAPOINTS WITHIN BIN, PER
        % RENDITION OF NONTARG
        numdatperrend = nrendsinbin/sum(indsthis); % number of datapoints per locked rend
        
        
        % -------------- 2) MIN NUM SONG BOUTS (based on unique time bins)
        if (1)
            tsongs = unique(tdev);
            nsongs = sum(tsongs>mintime_fromref & tsongs<maxtime_fromref);
            if nsongs<minsongs_inbin
                disp(['too few songs! sylcounter: ' num2str(ss)]);
                continue
            end
        end
        
        % ------------ 2) significant learning?
        siglearn = unique(DATBYREND.SigLearn(indsthis));
        assert(length(siglearn)<2);
        if onlyifSigLearn==1
            if siglearn==0
                continue
            end
        end
        
        % ========================= OUTPUT
        NtargCell{ss} = n_targ;
        NnontargCell{ss} = n_nontarg;
        FFbinnedCell{ss} = ffbinned_mean;
        TbinnedCell{ss} = tbinned;
        FFsinglebinMat(ss) = ffdev_singlebin;
        NumDatPerRend(ss) = numdatperrend;
        BirdNumMat(ss) = birdthis;
end



% %% ================== [COLLECT] RATIO OF HI TO LO DENS RENDS (TARG)
% N = cellfun(@mean, NtargCell);
% Nratio_hilo_targ = N(:,2)./N(:,1);
% 

Nratio_hilo_targ = nan;

