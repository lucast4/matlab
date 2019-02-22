function [OUTSTRUCT, OUTSTRUCT_XCOV, indsXcov] = lt_neural_POPLEARN_MatchOutstructs(OUTSTRUCT, OUTSTRUCT_XCOV, ...
    SwitchStruct, onlygoodexpt)

%% lt 2/20/19 - donwsamples OUTSTRUCT so that it matches OUTSTRUCT_XCOV
% i.e. it wil only contain cases (bird, expt, sw, motif) that also exist in
% OUTSTRUCT_XCOV.

% indsXcov is the corresponding ind in OUTSTRUCT_XCOV. (i.e. if multuple
% neuron pairs then saves all.

%% figure out the exact set of syls and expts to plot

if onlygoodexpt==1
    
    [OUTSTRUCT_XCOV] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
        SwitchStruct, 'xcov_spikes');
    
end

%% GO THRU OUTSTRUCT - ONLY KEEP CASES THAT HAVE EXACT MATCH IN OUTSTRUCT_XCOV
% need to do this becuase OUTSTRUCT contains pitch values..


indstokeep = [];
indsXcov = {};
indsXcov_all = {};
for i=1:length(OUTSTRUCT.bnum)
    
    % ----- look for corresponding ind in outstruc_xcov
    bnum = OUTSTRUCT.bnum(i);
    enum = OUTSTRUCT.enum(i);
    sw = OUTSTRUCT.switch(i);
    mnum = OUTSTRUCT.motifnum(i);
    
    indsthis = OUTSTRUCT_XCOV.bnum==bnum & OUTSTRUCT_XCOV.enum==enum & ...
        OUTSTRUCT_XCOV.switch==sw & OUTSTRUCT_XCOV.motifnum==mnum;
    
    if any(indsthis)
        indstokeep = [indstokeep; i];
        
        % ==== save the inds in xcov
        indsXcov = [indsXcov; find(indsthis)];
        
    end
    
    if any(indsthis)
    indsXcov_all = [indsXcov_all; find(indsthis)];
    else
        indsXcov_all = [indsXcov_all; {[]}];
    end
end

OUTSTRUCT.indsXcov_all = indsXcov_all;
% ===================== EXTRACT INDS THAT WANT TO KEEP
OUTSTRUCT = lt_structure_RmvEmptyField(OUTSTRUCT);
OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);