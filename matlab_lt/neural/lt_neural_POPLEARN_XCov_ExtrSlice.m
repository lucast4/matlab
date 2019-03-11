function [OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XCov_ExtrSlice(OUTSTRUCT_XCOV, ...
    OUTSTRUCT, PARAMS, twindow, alignto)
%% lt 2/1/19 - extracts scalars from desired t/f window
% CAN choose to be window relative to WN onset or to syl onset.

% ===== inds for time
t = PARAMS.xcenters_gram;
indst = t>=twindow(1) & t<=twindow(2);
assert(length(t)==size(OUTSTRUCT_XCOV.XcovgramBase{1},1), 'times are not matched to the data...');

%% Go thru each case, and extract
tmpBase = nan(length(OUTSTRUCT_XCOV.bnum), size(OUTSTRUCT_XCOV.XcovgramBase{1}, 2));
tmpWn = nan(length(OUTSTRUCT_XCOV.bnum), size(OUTSTRUCT_XCOV.XcovgramBase{1}, 2));

covslice_epochs = cell(length(OUTSTRUCT_XCOV.bnum), 1);
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
        assert(min(t)<=twindtmp(1) & max(t)>=twindtmp(2), 'not enough data to do this...');
    elseif strcmp(alignto, 'sylonset')
        % doi nothing, is default
    end
    
    % === collect slice
    tmpBase(i,:) = mean(covgram_base(indst, :), 1);
    tmpWn(i,:) = mean(covgram_wn(indst, :), 1);
    
    
    % ################# DO FOR XCOVGRA WITH MULTIPLE EPOCHS
%     covslice_epochs = nan(size(OUTSTRUCT_XCOV.XcovgramWN_epochs{1}, 3), ...
%         size(OUTSTRUCT_XCOV.XcovgramWN_epochs{1}, 2));
    if isfield(OUTSTRUCT_XCOV, 'XcovgramWN_epochs')
        xgram_epochs = OUTSTRUCT_XCOV.XcovgramWN_epochs{i};
        
        % ======== for each epoch, extract mean slice
        tmpepochs = squeeze(mean(xgram_epochs(indst, :,:),1))';
        covslice_epochs{i} = tmpepochs;
    end
end

% ======= save old versions
OUTSTRUCT_XCOV.XcovBase_orig = OUTSTRUCT_XCOV.XcovBase;
OUTSTRUCT_XCOV.XcovWN_orig = OUTSTRUCT_XCOV.XcovWN;

% ==== save new versions
OUTSTRUCT_XCOV.XcovBase = tmpBase;
OUTSTRUCT_XCOV.XcovWN = tmpWn;

OUTSTRUCT_XCOV.Xcovslice_epochs = covslice_epochs;

%% ======= save to params
PARAMS.Xcov_Slice_twindow = twindow;
PARAMS.Xcov_Slice_alignto = alignto;

