function [SwitchXCovStruct, PARAMS] = lt_neural_POPLEARN_XcovCalc(BirdsToSkip, ...
    binsize_spk, xcov_dattotake, xcovwindmax, normmethod, normspiketoprob, ...
    bregionpairtoget, removeIfLFPOnly, getXgram, windsize, windshift, PARAMS, ...
    SwitchCohStruct, SwitchStruct, MOTIFSTATS_pop, SummaryStruct, goodexptonly, ...
    dojitter, jitterwindSize)
%% ============= lt 2/2019 - does xcov analyses.

% BirdsToSkip = {'wh72pk12'};
% 
% % ====== binning spikes:
% binsize_spk = 0.0025; % default, 5ms bins for cross corr
% % xcov_dattotake = [-0.12 0.02]; % rel syl onset.
% % xcov_dattotake = [-0.05 0.02]; % rel syl onset.
% xcov_dattotake = [-0.1 0]; % rel syl onset. [BEST WINDOW]
% % xcov_dattotake = [-0.08 0.02]; % rel syl onset.
% % xcovwindmax = 0.06; % seconds
% xcovwindmax = 0.06; % seconds
% normmethod = 'unbiased';
% % normmethod = 'coeff';
% 
% % ======= norm bined spikes to prob in bin? (i.e. sum to 1)
% normspiketoprob = 0; % if 0, then uses spike count
% 
% % ======= bregionpairtoget
% bregionpairtoget = 'LMAN-RA';
% 
% % ======== to remove units that were extracted only for LFP
% removeIfLFPOnly = 1;
% 
% % ====== FOR GETTING RUNNING XCOV
% getXgram=1;
% windsize = 0.1;
% windshift = 0.005;
% % windsize = 0.08;
% % windshift = 0.005;
% % ---

%%

if goodexptonly ==1
    birdtoget = [1 2 4];
    swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
    firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
    % for a given switch did not have a previous switch on the same day
    % swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
    % firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
    % for a given switch did not have a previous switch on the same day
    firstswitchofday=1;
    indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
        firstswitchfortarget_withinday, firstswitchofday, birdtoget);
else
    indtoget_b_e_s  = [];
end
%% =============== PARAMS
windlist = [-0.13:windshift:(0.055-windsize)];
windlist = [windlist' windlist'+windsize];

% ====================== SAVE TO PARAMS
PARAMS.Xcov.xcov_dattotake = xcov_dattotake;
PARAMS.Xcov.binsize_spk = binsize_spk;
PARAMS.Xcov.normmethod = normmethod;
PARAMS.Xcov.bregionpairtoget = bregionpairtoget;
PARAMS.Xcov.Xgram.getXgram = getXgram;
PARAMS.Xcov.Xgram.windlist = windlist;


SwitchXCovStruct = struct;


%% ================== RUN
for i=1:length(SwitchCohStruct.bird)
    bname = SwitchStruct.bird(i).birdname;
    if any(ismember(BirdsToSkip, bname))
        continue
    end
    for ii=1:length(SwitchCohStruct.bird(i).exptnum)
        disp([i ii]);
        for iii=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
            DAT = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii);
            
            if ~isempty(indtoget_b_e_s)
                if ~ismember([i ii iii], indtoget_b_e_s, 'rows')
                    continue
                end
            end
            
            % ========= FOR THIS DATASET, GET SPIKES OVER ALL TRIALS FOR
            % RELEVANT CHANNELS.
            nmotifs = length(DAT.motifnum);
            for mm=1:nmotifs
                neurset = DAT.motifnum(mm).neursetused;
                
                if isempty(neurset)
                    continue
                end
                
                % ======== get segmentsextract data
                neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{neurset};
                bregionlist = {SummaryStruct.birds(i).neurons(neurlist).NOTE_Location};
                segextract_all = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(neurset).motif(mm).SegExtr_neurfakeID;
                segcommon = segextract_all(1).SegmentsExtract;
                chanlist = [SummaryStruct.birds(i).neurons(neurlist).channel];
                
                if removeIfLFPOnly==1
                    % -- check unit type for each neuron
                    isLFP = strcmp({SummaryStruct.birds(i).neurons(neurlist).NOTE_PutativeCellType}, 'LF');
                    neurlist(isLFP) = [];
                    bregionlist(isLFP) = [];
                    segextract_all(isLFP) = [];
                end
                
                if length(segcommon)<3 % too few trials
                    continue
                end
                
                % ========= GO THRU ALL PAIRS OF NEURONS. ONLY GET PAIRWISE
                % STATS IF THEY ARE BREGIONS OF INTEREST.
                ccRealAllPair = {};
                ccShiftAllPair = {};
                
                ccRealAuto = {};
                ccShiftAuto = {};
                
                
                ccRealAllPair_allwind= {};
                ccShiftAllPair_allwind = {};

                ccRealAllPair_Auto1_allwind = {};
                ccShiftAllPair_Auto1_allwind = {};
                ccRealAllPair_Auto2_allwind = {};
                ccShiftAllPair_Auto2_allwind = {};

                ccLags = [];
                neurPair = [];
                chanPair = [];
                % =========== TO COLLECT AU
                for nn=1:length(neurlist)
                    for nnn=nn+1:length(neurlist)
                        if strcmp([bregionlist{nn} '-' bregionlist{nnn}], bregionpairtoget)
                            seg1 = segextract_all(nn).SegmentsExtract;
                            seg2 = segextract_all(nnn).SegmentsExtract;
                            neur1 = neurlist(nn);
                            neur2 = neurlist(nnn);
                        elseif strcmp([bregionlist{nnn} '-' bregionlist{nn}], bregionpairtoget)
                            seg1 = segextract_all(nnn).SegmentsExtract;
                            seg2 = segextract_all(nn).SegmentsExtract;
                            neur1 = neurlist(nnn);
                            neur2 = neurlist(nn);
                        else
                            continue
                            % since is not desired pair
                        end
                        
                        
                        %% =============== COLLECT ALL PAIRWISE THINGS
                        % ======== 1) XCOV (SPIKES)
                        maxdur = min([segcommon.global_offtime_motifInclFlank] ...
                            - [segcommon.global_ontime_motifInclFlank]);
                        
                        seg1 = lt_neural_QUICK_SpkBinned(seg1, maxdur, ...
                            binsize_spk, 1, dojitter, jitterwindSize);
                        seg2 = lt_neural_QUICK_SpkBinned(seg2, maxdur, ...
                            binsize_spk, 1, dojitter, jitterwindSize);
                        
                        dattmp1 = struct;
                        dattmp1.SegmentsExtract = seg1;
                        dattmp2 = struct;
                        dattmp2.SegmentsExtract = seg2;
                        
                        [ccRealAll, ccShiftAll, lags_sec] = lt_neural_POP_GetXcov(...
                            dattmp1, dattmp2, xcov_dattotake, PARAMS.motif_predur, ...
                            xcovwindmax, binsize_spk, 0, 0, normmethod, normspiketoprob);
                        
                        ccRealAllPair = [ccRealAllPair; ccRealAll];
                        ccShiftAllPair = [ccShiftAllPair; ccShiftAll];
                        ccLags = lags_sec;
                        neurPair = [neurPair; [neur1 neur2]];
                        
                        chanPair = [chanPair; ...
                            [SummaryStruct.birds(i).neurons([neur1 neur2]).channel]];
                        
                        
                        %% ================== GET AUTOCOVARIANCE FUNCTIONS
                        % NOTE: DID NOT INSPECT THIS A SECOND TIME
                        % (2/18/20)
                        % AS WELL.
                        ccauto = {};
                        ccautoshift = {};
                        % -- 1) unit 1
                        [cc, ccshift] = lt_neural_POP_GetXcov(...
                            dattmp1, dattmp1, xcov_dattotake, PARAMS.motif_predur, ...
                            xcovwindmax, binsize_spk, 0, 0, normmethod, normspiketoprob);
                        ccauto = [ccauto cc];
                        ccautoshift = [ccautoshift ccshift];
                        
                        % -- 1) unit 2
                        [cc, ccshift] = lt_neural_POP_GetXcov(...
                            dattmp2, dattmp2, xcov_dattotake, PARAMS.motif_predur, ...
                            xcovwindmax, binsize_spk, 0, 0, normmethod, normspiketoprob);
                        ccauto = [ccauto cc];
                        ccautoshift = [ccautoshift ccshift];
                        
                        
                        ccRealAuto = [ccRealAuto; ccauto];
                        ccShiftAuto = [ccShiftAuto; ccautoshift];
                        
                        
                        %% ================ COLLECT XCOV, WITH SHIFTING TIME WINDOW
                        if getXgram==1
                            nwind = size(windlist,1);
                            ccreal_allwind = cell(1, nwind);
                            ccshift_allwind = cell(1, nwind);
                            for ww = 1:nwind
                                windthis = windlist(ww,:);
                                [ccreal, ccshift] = lt_neural_POP_GetXcov(...
                                    dattmp1, dattmp2, windthis, PARAMS.motif_predur, ...
                                    xcovwindmax, binsize_spk, 0, 0, normmethod, ...
                                    normspiketoprob);
                                
                                % -- save
                                ccreal_allwind{ww} = ccreal;
                                ccshift_allwind{ww} = ccshift;
                            end
                            % --- save output
                            ccRealAllPair_allwind = [ccRealAllPair_allwind; ccreal_allwind];
                            ccShiftAllPair_allwind = [ccShiftAllPair_allwind; ccshift_allwind];
                        end
                        
                        %% ================ COLLECT XCOV, WITH SHIFTING TIME WINDOW [AUTOCOVARIANCES]
                        if getXgram==1
                            nwind = size(windlist,1);
                            ccreal_allwind = cell(1, nwind);
                            ccshift_allwind = cell(1, nwind);
                            for ww = 1:nwind
                                windthis = windlist(ww,:);
                                [ccreal, ccshift] = lt_neural_POP_GetXcov(...
                                    dattmp1, dattmp1, windthis, PARAMS.motif_predur, ...
                                    xcovwindmax, binsize_spk, 0, 0, normmethod, ...
                                    normspiketoprob);
                                
                                % -- save
                                ccreal_allwind{ww} = ccreal;
                                ccshift_allwind{ww} = ccshift;
                            end
                            % --- save output
                            ccRealAllPair_Auto1_allwind = [ccRealAllPair_Auto1_allwind; ccreal_allwind];
                            ccShiftAllPair_Auto1_allwind = [ccShiftAllPair_Auto1_allwind; ccshift_allwind];
                        end
                        
                        %% ================ COLLECT XCOV, WITH SHIFTING TIME WINDOW [AUTOCOVARIANCES]
                        if getXgram==1
                            nwind = size(windlist,1);
                            ccreal_allwind = cell(1, nwind);
                            ccshift_allwind = cell(1, nwind);
                            for ww = 1:nwind
                                windthis = windlist(ww,:);
                                [ccreal, ccshift] = lt_neural_POP_GetXcov(...
                                    dattmp2, dattmp2, windthis, PARAMS.motif_predur, ...
                                    xcovwindmax, binsize_spk, 0, 0, normmethod, ...
                                    normspiketoprob);
                                
                                % -- save
                                ccreal_allwind{ww} = ccreal;
                                ccshift_allwind{ww} = ccshift;
                            end
                            % --- save output
                            ccRealAllPair_Auto2_allwind = [ccRealAllPair_Auto2_allwind; ccreal_allwind];
                            ccShiftAllPair_Auto2_allwind = [ccShiftAllPair_Auto2_allwind; ccshift_allwind];
                        end
                    end
                end
                
                
                % ================ OUTPUT FOR THIS MOTIF
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccRealAuto = ccRealAuto;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccShiftAuto = ccShiftAuto;

                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccRealAllPair_Auto1_allwind = ccRealAllPair_Auto1_allwind;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccShiftAllPair_Auto1_allwind = ccShiftAllPair_Auto1_allwind;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccRealAllPair_Auto2_allwind = ccRealAllPair_Auto2_allwind;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccShiftAllPair_Auto2_allwind = ccShiftAllPair_Auto2_allwind;

                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccRealAllPair = ccRealAllPair;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccShiftAllPair = ccShiftAllPair;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccRealAllPair_allwind = ccRealAllPair_allwind;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccShiftAllPair_allwind = ccShiftAllPair_allwind;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).ccLags = ccLags;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).neurPair = neurPair;
                SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm).chanPair = chanPair;
                
            end
        end
    end
end


%% ===== UPDATE PARAMS
PARAMS.Xcov_ccLags = lags_sec;    
