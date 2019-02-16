function [OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XCov_ExtrSlice(OUTSTRUCT_XCOV, ...
    OUTSTRUCT, PARAMS, twindow, alignto)
%% lt 2/1/19 - extracts scalars from desired t/f window
% CAN choose to be window relative to WN onset or to syl onset.

% ===== inds for time
t = PARAMS.xcenters_gram;
indst = t>=twindow(1) & t<=twindow(2);

%% Go thru each case, and extract
tmpBase = nan(size(OUTSTRUCT_XCOV.XcovBase));
tmpWn = nan(size(OUTSTRUCT_XCOV.XcovBase));
for i=1:length(OUTSTRUCT_XCOV.bnum)
    
    covgram_base = OUTSTRUCT_XCOV.XcovgramBase{i};
    covgram_wn = OUTSTRUCT_XCOV.XcovgramWN{i};
    
    if strcmp(alignto, 'wnonset')
        % --- recalcualte WN onset for this case
        % ============= FIGURE OUT WN TIME FOR THIS SWITCH
        bnum = OUTSTRUCT_XCOV.bnum(i);
        enum = OUTSTRUCT_XCOV.enum(i);
        sw = OUTSTRUCT_XCOV.switch(i);
        
        indstmp = OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & OUTSTRUCT.switch==sw;
        wntime = unique(OUTSTRUCT.WNonset(indstmp));
        assert(length(wntime)==1);
        
        twindtmp = twindow+wntime;
        indst = t>=twindtmp(1) & t<=twindtmp(2);
    else
        sadfasdf;
    end
    
    % === collect slice
    tmpBase(i,:) = mean(covgram_base(indst, :), 1);
    tmpWn(i,:) = mean(covgram_wn(indst, :), 1);
    
end

% ======= save old versions
OUTSTRUCT_XCOV.XcovBase_orig = OUTSTRUCT_XCOV.XcovBase;
OUTSTRUCT_XCOV.XcovWN_orig = OUTSTRUCT_XCOV.XcovWN;

% ==== save new versions
OUTSTRUCT_XCOV.XcovBase = tmpBase;
OUTSTRUCT_XCOV.XcovWN = tmpWn;

%% ======= save to params
PARAMS.Xcov_Slice_twindow = twindow;
PARAMS.Xcov_Slice_alignto = alignto;

