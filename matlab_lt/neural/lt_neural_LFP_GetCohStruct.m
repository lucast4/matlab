function [COHSTRUCT] = lt_neural_LFP_GetCohStruct(LFPSTRUCT, PARAMS, SummaryStruct, ...
    movingwin)
%% extracts and saves coherence for all trials.

% motif_predur = PARAMS.motif_predur;
% motif_postdur = PARAMS.motif_postdur;
PrePostRelSameTime = PARAMS.alignbyonset;
assert(~isempty(PrePostRelSameTime), 'then I have to check what empty means by default...');

% BirdsToPlot = {'pu69wh78', 'wh44wh39'};
% % SetsToSkip = {'1-2-2'};
% SetsToSkip = {};

%% COHERNECE PARAMS
lt_switch_chronux(1);

% movingwin = [0.1 0.01];

params = struct;
params.fpass = [1/movingwin(1) 200];
w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
% w = 20; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
tw = movingwin(1)*w;
params.tapers = [tw 2*tw-1];
params.Fs = 1500; % hard coded fs for LFP;


% lfpband_lo = 15;
% lfpband_hi = 35;
% 
% ffhilim = 90; % for plots;


% ======== for power spectra
% twind_spectrum = [-0.08 0]; % relative to syl onset

%% ===== make dir to save
savedir = ['/bluejay5/lucas/analyses/neural/LFP/PROCESSED/' PARAMS.savemarker];
mkdir(savedir);

%%
COHSTRUCT = struct;
% ==============
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    birdname = SummaryStruct.birds(i).birdname;
    disp(birdname);
    
    numexpts = length(LFPSTRUCT.bird(i).experiment);
    %     numexpts = length(MOTIFSTATS_pop.birds(i).exptnum);
    for ee = 1:numexpts
        disp(['expt ' num2str(ee)]);
        
        %         numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ee).Sets_neurons);
        numsets = length(LFPSTRUCT.bird(i).experiment(ee).setnum);
        
        % =============== FOR THIS BIRD, COLLECT ALL COHERENCE
        for ss = 1:numsets
            disp(['neuron set ' num2str(ss)]);
            
            dat = LFPSTRUCT.bird(i).experiment(ee).setnum(ss);
            if isempty(dat.motif)
                continue
            end
            nummotifs = length(dat.motif);
            
            
            % =============== FOR EACH MOTIF, COLLECT ALL TRIALS FOR ALL CHAN PAIRS
            for mm=1:nummotifs
                disp(['motif ' num2str(mm)]);
                
                
                % =============== EXTRACT COHERENCE FOR THIS MOTIF ACROSS ALL CHANS
                % AND TRIALS
                LFPall = dat.motif(mm).LFP_chanbytrial;
                Chanlist = dat.motif(mm).Chanlist;
                ntrials = size(LFPall, 2);
                if isempty(LFPall)
                    continue
                end
                
                motifpredur_plus_flank = -dat.motif(mm).t_relons(1) + ...
                    (dat.motif(mm).t_relons(2)-dat.motif(mm).t_relons(1)); % predur + flanking that was added when extract LFP
                
                rowcount = 1;
                CohAllTrials = cell(nchoosek(length(Chanlist),2), ntrials); % trial, then chan pair
                PhiAllTrials = cell(nchoosek(length(Chanlist),2), ntrials); % trial, then chan pair
                S12AllTrials = cell(nchoosek(length(Chanlist),2), ntrials); % trial, then chan pair
                S1AllTrials = cell(nchoosek(length(Chanlist),2), ntrials); % trial, then chan pair
                S2AllTrials = cell(nchoosek(length(Chanlist),2), ntrials); % trial, then chan pair
                
                Chanpairs = [];
                for c = 1:length(Chanlist)
                    for cc = c+1:length(Chanlist)
                        chan1 = Chanlist(c);
                        chan2 = Chanlist(cc);
                        Chanpairs = [Chanpairs; [chan1 chan2]];
                        % =========== go thru all trials getting coherence
                        dat1 = cell2mat(LFPall(c,:));
                        dat2 = cell2mat(LFPall(cc,:));
                        
                        % ===========
                        %                 [C,phi,S12,S1,S2,t,f] = cohgramc(dat1, dat2, movingwin, params);
                        [C,phi,S12,S1,S2,t,ffbins] = cohgramc(dat1, dat2, movingwin, params);
                        C = single(C);
                        phi = single(phi);
                        S12 = single(S12);
                        S1 = single(S1);
                        S2 = single(S2);                      
                        
                        % ======== SAVE INTO CELL ARRAY
                        tmp = squeeze(mat2cell(C, size(C,1), size(C,2), ones(1, size(C,3))))';
                        CohAllTrials(rowcount, :) = tmp;
                        
                        tmp = squeeze(mat2cell(phi, size(phi,1), size(phi,2), ones(1, size(phi,3))))';
                        PhiAllTrials(rowcount, :) = tmp;

                        tmp = squeeze(mat2cell(S12, size(S12,1), size(S12,2), ones(1, size(S12,3))))';
                        S12AllTrials(rowcount, :) = tmp;

                        tmp = squeeze(mat2cell(S1, size(S1,1), size(S1,2), ones(1, size(S1,3))))';
                        S1AllTrials(rowcount, :) = tmp;

                        tmp = squeeze(mat2cell(S2, size(S2,1), size(S2,2), ones(1, size(S2,3))))';
                        S2AllTrials(rowcount, :) = tmp;

                        rowcount = rowcount+1;
                    end
                end
                t_relons = t-motifpredur_plus_flank;
                
                
                if (0) % OLD VERSION - I verfieid that outputs are identical to the current version
                    % old version is 3x time (slower)
                   tic
                   segextract = MOTIFSTATS_pop.birds(i).exptnum(ee).DAT.setnum(ss).motif(mm).SegExtr_neurfakeID.SegmentsExtract;
                    neurthis = MOTIFSTATS_pop.birds(i).exptnum(ee).DAT.setnum(ss).motif(mm).SegExtr_neurfakeID.neurID_orig;
                    [CohAllTrials, Chanpairs, t_relons, ffbins] = lt_neural_Coher_GetAllMotifsChans(...
                        SummaryStruct, i, neurthis, segextract, Chanlist, motif_predur, motif_postdur,PrePostRelSameTime);
                    toc
                end
                
                % ============== COLLECT
                % filesize too large to keep in memory - save to disk and
                % note down path in this structure
                save([savedir '/Coh_bird' num2str(i) '_expt' num2str(ee) '_set' num2str(ss) '_mot' num2str(mm) '.mat'], ...
                    'CohAllTrials');
                save([savedir '/phi_bird' num2str(i) '_expt' num2str(ee) '_set' num2str(ss)  '_mot' num2str(mm) '.mat'], ...
                    'PhiAllTrials');
                save([savedir '/S12_bird' num2str(i) '_expt' num2str(ee) '_set' num2str(ss)  '_mot' num2str(mm) '.mat'], ...
                    'S12AllTrials');
                save([savedir '/S1_bird' num2str(i) '_expt' num2str(ee) '_set' num2str(ss)  '_mot' num2str(mm) '.mat'], ...
                    'S1AllTrials');
                save([savedir '/S2_bird' num2str(i) '_expt' num2str(ee) '_set' num2str(ss)  '_mot' num2str(mm) '.mat'], ...
                    'S2AllTrials');
                
                % ------- save overall structure
                COHSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).Chanpairs = Chanpairs;
                COHSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).t_relons = t_relons;
                COHSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).ffbins = ffbins;
                
            end
        end
    end
end

save([savedir '/COHSTRUCT.mat'], 'COHSTRUCT');
save([savedir '/movingwin.mat'], 'movingwin');
save([savedir '/params.mat'], 'params');
lt_switch_chronux(0);
