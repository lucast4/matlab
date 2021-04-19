classdef Neural < handle
    
    properties
        OUTSTRUCT_XCOV
        N
        SwitchStruct
        OUTSTRUCT
        PARAMS
        SwitchCohStruct
        MOTIFSTATS_pop
    end
    
    methods
        
        function obj = Neural()
            
        end
        
        function load_data(obj)
            % Automatic pipeline to run all loading and preprocessing steps
            
            % ==== LOAD
            % ===================== [PUBLICATION]
            % 1) all data, 2/5/19, up to gr48, RALMANLearn6
            load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_05Feb2019_2142.mat');
            load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_05Feb2019_2142.mat');
            load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_05Feb2019_2142.mat');
            load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/05Feb2019_2142/COHSTRUCT.mat');
            
            settmp = 1;
            if isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif)
                settmp=2;
                assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif));
            end
            PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).t_relons;
            PARAMS.ffbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).ffbins;
            PARAMS.ffbinsedges = [20 35 80 130]; % edges, to plot timecourse in frequency bands
            
            lt_neural_Coher_PreProcess;
            
            %%%%%%%%%%%%%%%%%
            savedir = ['/bluejay5/lucas/analyses/neural/POPLEARN/OUTSTRUCT_XCOV/' PARAMS.savemarker];
            
            % ===========  [PUBLICATION]
            savemarker = '60mswind_031919_3bins_hiloFF';  % latest, and 6 epochs for Xcovgram
            load([savedir '/SwitchXCovStruct_' savemarker '.mat']);
            
            % ====================== [SAME, but 4 bins for both xcovgram and epochs]
            savemarker = '60mswind_031919_4binsAll';  % 4 bins for both xcovgram and epochs
            load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV');
            load([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS')
            
            %%%%%%%%%%%%%%%% AD HOC COMBINING OF EPOCHS
            % == PU69 RL2, since there is gap in trials, need to combine.
            numepochs = size(OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{1}, 3);
%             disp(['nume pochs = ' num2str(numepochs) '?']);
%             pause;
            assert(size(OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs,2)==2, 'then need to modify code for OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs below');
            
            btoget = find(strcmp({SummaryStruct.birds.birdname}, 'pu69wh78'));
            etoget = find(strcmp({SwitchStruct.bird(btoget).exptnum.exptname}, 'RALMANlearn2'));
            
            indsthis = find(OUTSTRUCT_XCOV.bnum==btoget & OUTSTRUCT_XCOV.enum==etoget);
            
            for j=indsthis'
                
                % ====================
                tmp = OUTSTRUCT_XCOV.XcovgramWN_epochs{j};
                tmp(:,:, numepochs) = nanmean(tmp, 3);
                tmp(:,:, 1:numepochs-1) = nan;
                OUTSTRUCT_XCOV.XcovgramWN_epochs{j} = tmp;
                
                % ====================
                tmp = OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{j, 1};
                tmp(:,:, numepochs) = nanmean(tmp, 3);
                tmp(:,:, 1:numepochs-1) = nan;
                OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{j,1} = tmp;
                
                tmp = OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{j,2};
                tmp(:,:, numepochs) = nanmean(tmp, 3);
                tmp(:,:, 1:numepochs-1) = nan;
                OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{j,2} = tmp;
                
                % ====================
                tmp = OUTSTRUCT_XCOV.trialedges_epoch{j};
                tmp([numepochs numepochs+1]) = [tmp(1) tmp(end)];
                tmp(1:numepochs-1) = nan;
                OUTSTRUCT_XCOV.trialedges_epoch{j} = tmp;
            end
            
            %%%%%%%%%%%%%%%%%%
            % ========= GOOD: [60ms window]
            alignto = 'sylonset';
            twindows = {[-0.06 -0.03]}; % one array for each window [x centers...] [inclusive]
            lagwindows = {[-0.009 0.006]}; % OK
            OUTSTRUCT_XCOV = lt_neural_POPLEARN_XCov_ExtrScal(OUTSTRUCT_XCOV, OUTSTRUCT, ...
                PARAMS, twindows, lagwindows, alignto);

            % ======= GOOD [60ms window]
            twindow = [-0.06 -0.03]; % one array for each window [x centers...] [inclusive]
            alignto = 'sylonset';
            OUTSTRUCT_XCOV = lt_neural_POPLEARN_XCov_ExtrSlice(OUTSTRUCT_XCOV, ...
                OUTSTRUCT, PARAMS, twindow, alignto);

            %%%%%%%%%% PUT INTO SELF
            obj.initialize_data(OUTSTRUCT_XCOV, SwitchStruct, only_good_expt, ...
                OUTSTRUCT, PARAMS, SwitchCohStruct, MOTIFSTATS_pop)
        end
        
        function initialize_data(obj, OUTSTRUCT_XCOV, SwitchStruct, only_good_expt, ...
                OUTSTRUCT, PARAMS, SwitchCohStruct, MOTIFSTATS_pop)
            
            % Pass in already loaded and preprocessed data
            % ===== filter outstruct
            if only_good_expt
                [OUTSTRUCT_XCOV, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
                    SwitchStruct, 'xcov_spikes');
            end
            
            obj.OUTSTRUCT_XCOV = OUTSTRUCT_XCOV;
            obj.N = size(obj.OUTSTRUCT_XCOV.bnum,1);
            
            % outstruct,
            obj.OUTSTRUCT = OUTSTRUCT;
            
            obj.SwitchStruct = SwitchStruct;
            
            obj.PARAMS = PARAMS;
            obj.SwitchCohStruct = SwitchCohStruct;
            obj.MOTIFSTATS_pop = MOTIFSTATS_pop;
            
            % Preprocess
            obj.preprocess();
            
            % Sanity checks
            obj.check_preprocess();
            
        end
        
        function preprocess(obj)
            % =========== [EXTRACTION] GET LEARNING INTO BINS
            % === ALSO PLOTS SAMPLE SIZES...
            windowprem = [-0.1 0]; % for counting spikes
            obj.OUTSTRUCT_XCOV = lt_neural_POPLEARN_Xcov_ExtrLearn(obj.OUTSTRUCT, obj.OUTSTRUCT_XCOV, ...
                obj.SwitchStruct, obj.PARAMS, obj.SwitchCohStruct, obj.MOTIFSTATS_pop, windowprem);
            
            
        end
        
        function check_preprocess(obj)
            % Check sanity checks
            
            % 1) should all be LMAN-RA
            assert(all(strcmp(obj.OUTSTRUCT_XCOV.bregionpair, 'LMAN-RA')), 'assumes n1 is lMAN, n2 is RA');
            
            disp('Passed sanity checjks!');
        end
        
        %% General purpose functions
        
        function get_inds_grp(obj)
            
            
        end
        
        function cellarr = replace_nan_with_empty(obj, cellarr)
            % replaces all elements in cellarr that are nan with an empty
            % cell
            n = 0;
            for i=1:length(cellarr)
                if ~iscell(cellarr{i})
                    if isnan(cellarr{i})
                        cellarr{i} = [];
                        %                         n = n+1;
                    end
                end
            end
        end
        
        function dat = flip_if_adaptive(obj, dat)
            % Given dat, flips it based on direction of leanring.
            % INPUT:
            % - dat, N x 2, where first col is fflo trials, second is ffhi.
            
            learndirTarg = obj.OUTSTRUCT_XCOV.learndirTarg;
            dat(learndirTarg==-1, :) = fliplr(dat(learndirTarg==-1, :));
            
        end
        
        %% Extract data (single channel, n spikes)
        
        function [nspksOut, indshi, indslo, indskeep]= extract_single_channel_spikes(obj, ...
                area, epochthis)
            % Given area (LMAN, RA) and epoch (0=baseline, 1, 2, 3,..
            % for epochs, return spk rates, speartely for each trial
            % across all channels.
            % INPUT:
            % - no_repeat_chans, then returns smaller data, where each row
            % is a unique bird, expt, switch, channel, motif.
            % RETUNRS:
            % - nspksMat, cell (N,1), N num pair sof channels. each element K
            % x 1, where K is num trials. OR is empty, if no data for this
            % epoch.
            % - indshi, indslo, each cell array (N,1) of logical arrays (1,K),
            % where is 1 if that trial is hi or lo pitch. cell element will
            % be empty cell if no data for that row,
            % - original indicces in outstruct, to allow going back.
            % - leandirTarg, array, N,1, where -1 or 1 dep on learning
            % direction for target syl for that ept. no nans.
            % - indskeep, logical arrya, N,1, which can use to filter to
            % keep only unqiue unique bird, expt, switch, channel, motif.
            
            OUTSTRUCT_XCOV = obj.OUTSTRUCT_XCOV;
            
            epochthis = epochthis +1; % since 0 is baseline.
            ind_area = find(strcmp(area, {'LMAN', 'RA'}));
            
            % == old version, doesn't work if epochthis is array length>1
            nspksAll= cellfun(@(x)x(epochthis).nspksByNeuron, OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
            
            % extract for this region. use for loop to take care of any
            % nan, which due to this epoch not having any data.
            nspksOut = cell(size(nspksAll,1), 1);
            for i=1:length(nspksAll)
                if iscell(nspksAll{i})
                    nspksOut{i} = nspksAll{i}{ind_area};
                else
                    assert(isnan(nspksAll{i}));
                end
            end
            
            % count how many empty
            n = find(cellfun(@(x)isempty(x), nspksOut));
            disp('num rows empty (no dat for epoch');
            disp(length(n));
            
            % Also get the inds for hi and lo pitch, and learn dir, which
            % can be used to say if is adaptive of now.
            indshi = cellfun(@(x)[x(epochthis).inds_hi], OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
            indslo = cellfun(@(x)[x(epochthis).inds_lo], OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
            % convert nans to empty
            indshi = obj.replace_nan_with_empty(indshi);
            indslo = obj.replace_nan_with_empty(indslo);
            
            % find the unique channels
            [indsgrp, indsgrpU, X_cell] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
                OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch, ...
                OUTSTRUCT_XCOV.chanpair_actual(:,ind_area), OUTSTRUCT_XCOV.motifID_unique});
            % sanity check, for each unique grp, confirm that spk rates are
            % identical across all inds for that grp
            indskeep = zeros(length(indsgrp),1);
            
            tmptmp = [];
            for i=1:length(indsgrpU)
                grpthis = indsgrpU(i);
                indsthis = find(indsgrp==grpthis);
                tmptmp = [tmptmp; indsthis];
                matthis = cell2mat(cellfun(@(x)transpose(x), nspksOut(indsthis), 'UniformOutput', false)); % n x numtrials.
                if size(matthis,1)>1
                    tmp = diff(matthis);
                    %                 disp(tmp);
                    %                 disp(grpthis);
                    assert(all(tmp(:)==0));
                else
                    % ok, only one datapoint, must be unique
                end
                
                %                 % make sure no empty
                %                 if iscell(matthis)
                %                     disp(matthis);
                %                     assert(false)
                %                 end
                
                if isempty(matthis)
                    dum =cellfun(@(x)transpose(x), nspksOut(indsthis), 'UniformOutput', false);
                    assert(all(cellfun(@(x)isempty(x), dum)), 'some trials empty some not? why? and should make indskeep the one not empty');
                    %                     disp(cellfun(@(x)transpose(x), nspksOut(indsthis), 'UniformOutput', false));
                    %                     assert(false)
                end
                
                if any(isnan(matthis(:)))
                    disp(matthis)
                    assert(false, 'no nans should be present at all, since converted to empties')
                end
                
                % keep the first ind
                indskeep(indsthis(1))=1;
            end
            
            assert(length(tmptmp)==length(nspksOut), 'somehow skipped some inds...');
            
            disp('Num unique datapoints found:');
            disp(sum(indskeep));
        end
        
        
        function [nspksOut, indshi, indslo, indskeep] = extract_single_channel_spikes_mean(obj, area, epochlist)
            % Same as extract_single_channel_spikes, but concatenates over
            % epochs. is ok if any epoch is nan, just skips.
            % INPUT:
            % - epochlist, 0 (baseline), 1, 2, 3,...
            
            dat = struct;
            for i=1:length(epochlist)
                epoch = epochlist(i);
                [nspksOut, indshi, indslo, indskeep] = obj.extract_single_channel_spikes(area, epoch);
                dat(i).nspks = nspksOut;
                dat(i).indshi = indshi;
                dat(i).indslo = indslo;
                dat(i).indskeep = indskeep;
            end
            
            % convert all to doubple
            for i=1:length(nspksOut)
                %                disp(nspksOut{i});
                nspksOut{i} = double(nspksOut{i});
            end
            
            % check that all indskeep are identcial
            for i=2:length(dat)
                assert(all(dat(i).indskeep==dat(i-1).indskeep));
            end
            indskeep = dat(1).indskeep;
            
            % concatenate
            nspksOut = cell(length(dat(1).nspks), 1);
            indshi = cell(length(dat(1).nspks), 1);
            indslo = cell(length(dat(1).nspks), 1);
            
            for i=1:length(dat)
                for j=1:length(nspksOut)
                    
                    nspksOut{j} = [nspksOut{j}; dat(i).nspks{j}];
                    indshi{j} = [indshi{j}, dat(i).indshi{j}];
                    indslo{j} = [indslo{j}, dat(i).indslo{j}];
                end
            end
        end
        
        function [dat, indskeep] = extract_dat_channel_spikes(obj, area, ...
                epochlist, splitver)
            % Extract clean output,
            % INPUT:
            % - splitver, {'adaptive', 'ff'}, wiether to make cols be
            % nonadapt, adapt or lo hi.
            % dat array, (N,2), where N is total rows, and cols are [lo,
            % hi]
            % indskeep, inds to use if want unique chans
            
            [nspksOut, indshi, indslo, indskeep] = ...
                obj.extract_single_channel_spikes_mean(area, epochlist);
            
            % take mean over inds
            dat = nan(length(nspksOut), 2);
            %             tmp = cellfun(@(x,y)x(logical(y)), nspksOut, indslo, 'UniformOutput', false);
            %             dat(:,1) = cellfun(@(x)nanmean(x), tmp, 'UniformOutput', true);
            %             disp(nspksOut{25})
            %             disp(nspksOut{26})
            %             assert false
            dat(:,1) = cellfun(@(x,y)double(nanmean(x(logical(y)))), nspksOut, indslo, 'UniformOutput', true);
            %             tmp = cellfun(@(x,y)x(logical(y)), nspksOut, indshi, 'UniformOutput', false);
            %             for i=1:length(tmp)
            %                 disp(i)
            %                 disp(tmp{i});
            %                 nanmean(tmp{i});
            %             end
            dat(:,2) = cellfun(@(x,y)double(nanmean(x(logical(y)))), nspksOut, indshi, 'UniformOutput', true);
            
            if strcmp(splitver, 'adaptive')
                dat = obj.flip_if_adaptive(dat);
                %                 learndirTarg = obj.OUTSTRUCT_XCOV.learndirTarg;
                %                 dat(learndirTarg==-1, :) = fliplr(dat(learndirTarg==-1, :));
            else
                assert(strcmp(splitver, 'ff'));
            end
            
        end
        
        %% Extract - xcov
        
        function dat = extract_xcov_ffslits_singleepoch(obj, epoch)
            % epoch 0: base
            % epoch 1, 2, ... WN
            % RETURNS:
            % - dat, N x 2, where cols are [fflo, ffhi]
            
            OUTSTRUCT_XCOV = obj.OUTSTRUCT_XCOV;
            
            if epoch==0
                dat = OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit;
                dat = cellfun(@(x)transpose(x), dat, 'UniformOutput', false);
                dat = cell2mat(dat);
            else
                dat = cellfun(@(x)x{1}(epoch,:), ...
                    OUTSTRUCT_XCOV.xcovscalEpochs_window_FFsplit, 'UniformOutput', false);
                dat = cell2mat(dat);
            end
        end
        
        function [dat] = extract_dat_xcov_ffsplits(obj, epochlist, splitver)
            % Combines across multiple epochs (nanmean)
            % INPUTS:
            % - epochlist, where 0 is baseline, 1, 2, 3, ... are wn
            % - splitver, {'adaptive', 'ff'}
            % --- Can choose if reorder based on dir learning [nonad, ad],
            % or [fflo, ffhi]
            % RETURNS:
            % - dat, N x 2.
            
            dat = nan(obj.N, 2, length(epochlist)); % (N, 2, length(epochlist))
            for i=1:length(epochlist)
                epoch = epochlist(i);
                datthis = obj.extract_xcov_ffslits_singleepoch(epoch);
                dat(:,:, i) = datthis;
            end
            
            % take mean over epochs
            dat = nanmean(dat, 3);
            
            % flip any?
            if strcmp(splitver, 'adaptive')
                dat = obj.flip_if_adaptive(dat);
                %                 learndirTarg = obj.OUTSTRUCT_XCOV.learndirTarg;
                %                 dat(learndirTarg==-1, :) = fliplr(dat(learndirTarg==-1, :));
            else
                assert(strcmp(splitver, 'ff'));
            end
            
        end
        
        %% Extract other things
        function TrainDurs = plot_training_duration(obj)
            
            OUTSTRUCT = obj.OUTSTRUCT;
            OUTSTRUCT_XCOV = obj.OUTSTRUCT_XCOV;
            SwitchStruct = obj.SwitchStruct;
            onlygoodexpt = true;
            
            % ======= first downsample OUTSRRUCT to get only the experiments in xcov
            [OUTSTRUCT, OUTSTRUCT_XCOV, indsXcov] = lt_neural_POPLEARN_MatchOutstructs(...
                OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, onlygoodexpt);
            
            % ====== for each case compute z-scored laerning
            %             LearnTargDir_Z = nan(size(OUTSTRUCT.bnum));
            %             LearnTargDir_Slope = nan(size(OUTSTRUCT.bnum));
            %             TrainDuration = nan(size(OUTSTRUCT.bnum)); % time from mid of base bin to mid of training
            TrainDurs = nan(size(OUTSTRUCT.bnum));
            for i=1:length(OUTSTRUCT.bnum)
                
                ffvals = OUTSTRUCT.ffvals{i};
                tvals = OUTSTRUCT.tvals{i};
                
                if all(isnan(ffvals))
                    continue
                end
                
                % ---- USE THE INDS FROM XCOV
                %     indsbase = OUTSTRUCT_XCOV.inds_base_epoch{indsXcov{i}(1)};
                indsbase = OUTSTRUCT_XCOV.inds_base_allgood{indsXcov{i}(1)}; % take first, since all shuld be identical.
                indsWN = OUTSTRUCT_XCOV.inds_WN_epoch{indsXcov{i}(1)};
                
                % --- sanity check....
                assert(all(OUTSTRUCT.indsXcov_all{i} == indsXcov{i}));
                
                % get time of final trial
                traindur = (tvals(max(indsWN)) -  tvals(max(indsbase)))*24;
                TrainDurs(i) = traindur;
            end
            
            % get unique traindurs
            tmp = unique(TrainDurs);
            tmp = tmp(~isnan(tmp));
            TrainDurs = tmp;
            
            lt_figure(); hold on;
            plot(TrainDurs, 1, 'ok');
            
            xlabel(['train dur (last base rend to last WN rend)']);
            ylabel(['min: ' num2str(min(TrainDurs)), 'max: ' num2str(max(TrainDurs))]);
            title(['mean, std: ' num2str([mean(TrainDurs), std(TrainDurs)])]);
        end
        
        %% to extract summarized dat of different kinds.
        %         function dat = extract_dat_single_channel_spikes(obj, baseline_epochs, wn_epochs)
        %             % dat, struct, one row for each channel, with fields:
        %             % mean spk rate (mean rate over trials)
        %             % baseline, wn, adaptive, nonadaptive, LMAN, RA
        %
        %             % Collect data for baseline epochs.
        %             for i=1:length(arealist)
        %                 areathis = arealist{i};
        %
        %                 % baseline
        %                 epochlist = baseline_epochs;
        %                 [nspksOut, indshi, indslo, indskeep] = ...
        %                     obj.extract_single_channel_spikes_mean(areathis, epochlist);
        %
        %             end
        %
        %         end
        
        %% CHECK THINGS
        % e.g,, about params, things I might want to verify
        function check_what_frac_trials_used(obj)
            % For baseline and WN, what/which frac trials used?
            
            OUTSTRUCT_XCOV = obj.OUTSTRUCT_XCOV;
            
            i=1;
            inds_WN_epoch = OUTSTRUCT_XCOV.inds_WN_epoch{i};
            inds_WN_allgood = OUTSTRUCT_XCOV.inds_WN_allgood{i};
            
            disp('WN inds used ... WN inds all')
            disp(length(inds_WN_epoch));
            disp(length(inds_WN_allgood));
            
            inds_base_epoch = OUTSTRUCT_XCOV.inds_base_epoch{i};
            inds_base_allgood = OUTSTRUCT_XCOV.inds_base_allgood{i};
            
            disp('WN inds used ... WN inds all')
            disp(inds_base_epoch);
            disp(inds_base_allgood);
            
        end
        
        %% HELPERS
        function Y = reshape_by_syltype(obj, dat)
            % Given dat (N,K) reshapes into cell array {targ, same, diff}
            % where each array in cell is diff num rows (same num columns).
            
            istarg = obj.OUTSTRUCT_XCOV.istarg;
            issame = obj.OUTSTRUCT_XCOV.issame;
            
            Y = cell(1,3);
            Y{1} = dat(istarg==1, :);
            Y{2} = dat(istarg==0 & issame==1, :);
            Y{3} = dat(istarg==0 & issame==0, :);
        end
        
        %% PLOTS - general
        
        function plotDatFFsplits(obj, dat_base, dat_wn)
            % General purpose, plot data, with different stratifications
            % dat_base, dat_wn, are N,2 arrays, where columns are [fflo,
            % ffhi]
            
            % summary version, everything normalized
            dat_wnminbase = dat_wn - dat_base; % subtract base
            dat_wnminbase_himinlo = dat_wnminbase(:,2) - dat_wnminbase(:,1); % hi ff - lo ff
            Y = obj.reshape_by_syltype(dat_wnminbase_himinlo);
            obj.plotDistsUnpaired(Y);
        end
        
        function plotDistsUnpaired(obj, Y)
            % Y is cell, (1,M) where each element is a single column
            % vector.
            
            lt_figure; hold on;
            x = 1:size(Y,2);
            lt_plot_MultDist(Y, x, 1, 'k');
            
            % for each pair, compute p val
            for i=1:length(Y)
                for ii=i+1:length(Y)
                    x1 = Y{i};
                    x2 = Y{ii};
                    p = ranksum(x1, x2);
                    if p<0.5
                        lt_plot_pvalue(p, [num2str(i) 'vs' num2str(ii)], 1);
                    end
                end
            end
        end
        
        function plotDistsPaired(obj, dat)
            % dat is either:
            % --- array, N x M, where each row is assumed paired.
            % plots each col as s ingle distribution, and will plot
            % lines between dots across rows
            lt_figure; hold on;
            
            Y = {};
            for i=1:size(dat,2)
                Y{i} = dat(:,i);
            end
            x = 1:size(Y,2);
            lt_plot_MultDist(Y, x, 1, 'k');
            
            % plot lines
            if size(dat,2)>1
                for i=1:size(dat,1) % each row
                    plot(x, dat(i,:), '-ob');
                end
            end
            
            % if paired, then do paired test
            if size(dat,2)==2
                p = signrank(dat(:,1), dat(:,2));
                lt_plot_pvalue(p, 'signrank vs', 1);
            end
            
            n = size(dat,1);
            lt_plot_annotation(1, ['n=' num2str(n)], 'r');
        end
        
        
        %% PLOTS - specific
        function plot_xcov_ffsplits(obj)
            epochlist_base = [0];
            epochlist_wn = [3 4];
            splitver = 'adaptive';
            
            dat_base = obj.extract_dat_xcov_ffsplits(epochlist_base, splitver);
            dat_wn = obj.extract_dat_xcov_ffsplits(epochlist_wn, splitver);
            
            obj.plotDatFFsplits(dat_base, dat_wn);
            
        end
        
        
        function plot_single_channel_spikes(obj, baseline_epochs, ...
                wn_epochs)
            % Overview, summary plots of single-area spike rates.
            
            splitverlist = {'ff', 'adaptive'};
            arealist = {'LMAN', 'RA'};
            targsyl = 1;
            
            function plots(dat, indskeep, splitver, area)
                
                % 1) Plot all 4 datasets
                datthis = dat(indskeep==1);
                obj.plotDistsPaired(datthis);
                if strcmp(splitver, 'adaptive')
                    xlabel('BASE (non - ad) -- WN (non, ad)');
                elseif strcmp(splitver, 'ff')
                    xlabel('BASE (loff, hiff) -- WN (loff, hiff)');
                end
                ylabel('fr (hz)')
                title(['units x syl, ' area]);
                
                % 2) Plto norm to base
                % keep only unique chans
                datthis = dat(indskeep==1 & obj.OUTSTRUCT_XCOV.istarg==targsyl, :);
                % normalize to baseline
                datthis_norm = nan(size(datthis,1), 2);
                datthis_norm(:,1) = datthis(:,3) - datthis(:,1);
                datthis_norm(:,2) = datthis(:,4) - datthis(:,2);
                obj.plotDistsPaired(datthis_norm);
                if strcmp(splitver, 'adaptive')
                    xlabel('nonadaptive - adaptive (norm to base)');
                elseif strcmp(splitver, 'ff')
                    xlabel('lowff - hiff (norm to base)');
                end
                ylabel('fr (hz)')
                title(['targsyl ' area]);
            end
            
            for i=1:length(splitverlist)
                splitver = splitverlist{i};
                for j=1:length(arealist)
                    area = arealist{j};
                    % get baseline
                    [dat_base, indskeep_base] = obj.extract_dat_channel_spikes(area, ...
                        baseline_epochs, splitver);
                    % get WN
                    [dat_wn, indskeep_wn] = obj.extract_dat_channel_spikes(area, ...
                        wn_epochs, splitver);
                    assert(all(indskeep_base==indskeep_wn));
                    % Collect
                    indskeep = indskeep_base;
                    dat = [dat_base, dat_wn];
                    
                    plots(dat, indskeep, splitver, area);
                end
            end
        end
        
        %% load state
        
        
        
    end
end