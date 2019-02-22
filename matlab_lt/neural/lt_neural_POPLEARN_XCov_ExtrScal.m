function [OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XCov_ExtrScal(OUTSTRUCT_XCOV, OUTSTRUCT, ...
    PARAMS, twindows, lagwindows, alignto)
%% lt 2/1/19 - extracts scalars from desired t/f window
% CAN choose to be window relative to WN onset or to syl onset.

assert(length(PARAMS.xcenters_gram) == size(OUTSTRUCT_XCOV.XcovgramBase{1},1));
assert(length(PARAMS.Xcov_ccLags) == size(OUTSTRUCT_XCOV.XcovgramBase{1},2));

%% Go thru each case, and extract
tmp = cell(length(OUTSTRUCT_XCOV.bnum),1);
for i=1:length(OUTSTRUCT_XCOV.bnum)
    
    covgram_base = OUTSTRUCT_XCOV.XcovgramBase{i};
    covgram_wn = OUTSTRUCT_XCOV.XcovgramWN{i};
    
    % ============= FIGURE OUT WN TIME FOR THIS SWITCH
    bnum = OUTSTRUCT_XCOV.bnum(i);
    enum = OUTSTRUCT_XCOV.enum(i);
    sw = OUTSTRUCT_XCOV.switch(i);
    indstmp = OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & OUTSTRUCT.switch==sw;
    wntime = unique(OUTSTRUCT.WNonset(indstmp));
    assert(length(wntime)==1);
    
    
    % ===== go thru each window and collect base and WN
    xcovscal_window_epoch = nan(length(twindows), 2); % i.e. (num windows x 2 (base, WN));
    for ii=1:length(twindows)
        twind = twindows{ii};
        lagwind = lagwindows{ii};
        
        t = PARAMS.xcenters_gram;
        lags = PARAMS.Xcov_ccLags;
        
        % ===
        if strcmp(alignto, 'sylonset')
            indst = t>=twind(1) & t<=twind(2);
        elseif strcmp(alignto, 'wnonset')
            twindtmp = twind+wntime;
            assert(min(t)<twindtmp(1) & max(t)>twindtmp(2));
            indst = t>=twindtmp(1) & t<=twindtmp(2);
        else
            sadfasdf;
        end
        indslag = lags>=lagwind(1) & lags<=lagwind(2);
        
        % === collect scalar
        Y = nan(1,2); % base, WN
        Y(1) = mean(mean(covgram_base(indst, indslag)));
        Y(2) = mean(mean(covgram_wn(indst, indslag)));
        
        
        % ====== save output
        xcovscal_window_epoch(ii, :) = Y;
    end
    
    tmp{i} = xcovscal_window_epoch;
end

OUTSTRUCT_XCOV.Xcovscal_window_BaseWN = tmp;

%% ======= save to params
PARAMS.Xcov_Scalar_twindows = twindows;
PARAMS.Xcov_Scalar_lagwindows = lagwindows;
PARAMS.Xcov_Scalar_alignto = alignto;

