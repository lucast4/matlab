function [OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XCov_ExtrScal(OUTSTRUCT_XCOV, OUTSTRUCT, ...
    PARAMS, twindows, lagwindows, alignto)
%% lt 2/1/19 - extracts scalars from desired t/f window
% CAN choose to be window relative to WN onset or to syl onset.

assert(length(PARAMS.xcenters_gram) == size(OUTSTRUCT_XCOV.XcovgramBase{1},1));
assert(length(PARAMS.Xcov_ccLags) == size(OUTSTRUCT_XCOV.XcovgramBase{1},2));

%% Go thru each case, and extract
tmp = cell(length(OUTSTRUCT_XCOV.bnum),1);
tmp_epochs = cell(length(OUTSTRUCT_XCOV.bnum),1);
for i=1:length(OUTSTRUCT_XCOV.bnum)
    
    covgram_base = OUTSTRUCT_XCOV.XcovgramBase{i};
    covgram_wn = OUTSTRUCT_XCOV.XcovgramWN{i};
    
    % ============= FIGURE OUT WN TIME FOR THIS SWITCH
    bnum = OUTSTRUCT_XCOV.bnum(i);
    enum = OUTSTRUCT_XCOV.enum(i);
    sw = OUTSTRUCT_XCOV.switch(i);
    indstmp = OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & OUTSTRUCT.switch==sw;
    
    
    % ########################### EPOCHS
         if isfield(OUTSTRUCT_XCOV, 'XcovgramWN_epochs')
   xgram_epoch = OUTSTRUCT_XCOV.XcovgramWN_epochs{i};
   Nepoch = size(xgram_epoch,3);
        end
    
        
    % ===== go thru each window and collect base and WN
    xcovscal_window_epoch = nan(length(twindows), 2); % i.e. (num windows x 2 (base, WN));
         if isfield(OUTSTRUCT_XCOV, 'XcovgramWN_epochs')
    xcovscal_WNepochs = nan(length(twindows), Nepoch);
         end
    for ii=1:length(twindows)
        twind = twindows{ii};
        lagwind = lagwindows{ii};
        
        t = PARAMS.xcenters_gram;
        lags = PARAMS.Xcov_ccLags;
        
        % ===
        if strcmp(alignto, 'sylonset')
            indst = t>=twind(1) & t<=twind(2);
        elseif strcmp(alignto, 'wnonset')
            wntime = unique(OUTSTRUCT.WNonset(indstmp));
            assert(length(wntime)==1);
            
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

        
        
        % ==== colllect scalar for epochs
        if isfield(OUTSTRUCT_XCOV, 'XcovgramWN_epochs')
            
            xcovscal_epochs = squeeze(mean(mean(xgram_epoch(indst, indslag, :),1),2));
            xcovscal_WNepochs(ii, :) = xcovscal_epochs;
        end
    end
    
    tmp{i} = xcovscal_window_epoch;
    
    if isfield(OUTSTRUCT_XCOV, 'XcovgramWN_epochs') 
    tmp_epochs{i} = xcovscal_WNepochs;
    end
    
end

OUTSTRUCT_XCOV.Xcovscal_window_BaseWN = tmp;
if isfield(OUTSTRUCT_XCOV, 'XcovgramWN_epochs')
    OUTSTRUCT_XCOV.Xcovscal_window_Epochs = tmp_epochs;
end
%% ======= save to params
PARAMS.Xcov_Scalar_twindows = twindows;
PARAMS.Xcov_Scalar_lagwindows = lagwindows;
PARAMS.Xcov_Scalar_alignto = alignto;

