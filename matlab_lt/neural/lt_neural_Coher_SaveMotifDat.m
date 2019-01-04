function lt_neural_Coher_SaveMotifDat(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, ...
    PARAMS, save_xchan_means)
%% lt 1/3/19 - [SAVE SUMMARY, FOR COMPARISON ACROSS EXTRACTED DATASETS]
% === Saves once for each bird/expt/switch -
% === Therefore b/e/s must be matched across strtuctures that are comparing
% === saves: 1) mean coherogram (pre and post) and 2) scalar (pre
% and post), for each syllable (once for each channel pair). 3) is target?
% 4) same type? 5) motif name

% This is useful if, e.g., want to comapre the same motif/expt across
% different MOTIFSTATS (e.g. data vs. controls, or diff data extractions
% with different subsets of birds)... Useful becuase MOTIFSTATS And the
% other extraction structure are large memroy./

%% save dir
savedir = '/bluejay0/bluejay2/lucas/analyses/neural/COHERENCE/SaveMotifDat';

%% ============ get diff in coh struct for each case
for i=1:length(OUTSTRUCT.CohMean_Base)
    
    OUTSTRUCT.CohMean_WNminusBase{i} = OUTSTRUCT.CohMean_WN{i} - OUTSTRUCT.CohMean_Base{i};
    
end

%% =========== first take mean across channel pairs, if so desired.
if (0)
    fieldtoget = 'CohMean_WNminusBase';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
    
    
    
    fieldtoget = 'CohMean_WN';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat2] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
    
    fieldtoget = 'CohMean_Base';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat1] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
end
%% RUN

% --- unique switches
[indsgrp_switch, indsgrp_switch_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
% --- uqniue chan pairs
[indsgrp_chanpair, indsgrp_chanpair_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, ...
    OUTSTRUCT.chanpair});

for i=1:length(indsgrp_switch_unique)
    swgrpthis = indsgrp_switch_unique(i);
    bnum = unique(OUTSTRUCT.bnum(indsgrp_switch==swgrpthis));
    enum = unique(OUTSTRUCT.enum(indsgrp_switch==swgrpthis));
    swnum = unique(OUTSTRUCT.switch(indsgrp_switch==swgrpthis));
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    %     % ====== plot each channel pair its own line
    %     for chanpair = indsgrp_chanpair_unique'
    %         indsthis = indsgrp_switch==swgrpthis & indsgrp_chanpair==chanpair;
    %         if ~any(indsthis)
    %             continue
    %         end
    %         motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
    %         chnums = OUTSTRUCT.chanpair(indsthis,:);
    %         chnums = unique(chnums)'; assert(length(chnums)==2);
    %         cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    %         %         cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
    %
    %         assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
    %         % ---- sort in order of motifs
    %         [~, indsort] = sort(motifID);
    %         motifID = motifID(indsort);
    %         cohscal = cohscal(indsort);
    %         plot(motifID, cohscal, 'o-k');
    %         lt_plot_text(motifID(end)+0.3, cohscal(end), num2str(chnums), 'm', 9);
    %     end
    %     lt_plot_zeroline;
    %
    
    % ====== overall
    indsthis = indsgrp_switch==swgrpthis;
    
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    motifs = OUTSTRUCT.motifname(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
    %     cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
    cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis)');
    
    if save_xchan_means==1
        % --- grp coh matrix
        nrow = size(cohmat,1);
        ncol = size(cohmat,2);
        cohmatout = nan(size(cohmat,1), size(cohmat,2), length(unique(motifID)));
        
        for rr=1:nrow
            for cc=1:ncol
                cohmatout(rr, cc, :) = grpstats(squeeze(cohmat(rr, cc, :)), motifID);
            end
        end
        cohmat = cohmatout;
        
        % --- grp scalar
        cohscal = grpstats(cohscal, motifID, {'mean'});
        istarg = grpstats(istarg, motifID, {'mean'});
        issame = grpstats(issame, motifID, {'mean'});
        
        % --- unique motifs
        motifID = unique(motifID);
        
        %     lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
    else
        disp('HAVE NOT CODED FOR save_xchan_means==0')
        return;
    end
    
    
    % ##################### SAVE, FOR THIS SWITCH
    savestruct = struct;
    savestruct.cohmat_diff = cohmat;
    savestruct.cohscal_diff = cohscal;
    savestruct.istarg = istarg;
    savestruct.issame = issame;
    savestruct.motifID = motifID;
    
    fname = [savedir '/' bname '_' ename '_sw' num2str(swnum) '_' PARAMS.savemarker '.mat'];
    save(fname, 'savestruct');
    disp(['SAVED: ' fname]);
end
