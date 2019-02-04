function [OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XCov_ExtrScal(OUTSTRUCT_XCOV, OUTSTRUCT, ...
    PARAMS, twindows, lagwindows, alignto)
%% lt 2/1/19 - extracts scalars from desired t/f window 
% CAN choose to be window relative to WN onset or to syl onset.


%% Go thru each case, and extract 
tmp = cell(length(OUTSTRUCT_XCOV.bnum),1);
for i=1:length(OUTSTRUCT_XCOV.bnum)
    
    covgram_base = OUTSTRUCT_XCOV.XcovgramBase{i};
    covgram_wn = OUTSTRUCT_XCOV.XcovgramWN{i};
    
    % ===== go thru each window and collect base and WN
    xcovscal_window_epoch = nan(length(twindows), 2); % i.e. (num windows x 2 (base, WN));
    for ii=1:length(twindows)
       twind = twindows{ii};
       lagwind = lagwindows{ii};
       
       t = PARAMS.xcenters_gram;
       lags = PARAMS.Xcov_ccLags;
       
       % === 
       indst = t>=twind(1) & t<=twind(2);
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

