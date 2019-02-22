function [OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XcovSmth(...
    OUTSTRUCT_XCOV, PARAMS, dosmooth_sigma)
%% 2/18/19 - LT, smoothed final xcov (i.e. dat minus shuffle)
% after doing all xcov stuff, smooth by 1) interpolating, then smoothing
% Output binsize will always be set to 0.001 sec

binsize = PARAMS.Xcov_ccLags(2)-PARAMS.Xcov_ccLags(1);
if mod(binsize, 0.001)>0.0001
    binsize_final = 0.0005;
else
    binsize_final = 0.001;
end


xbins = PARAMS.Xcov_ccLags; % before interpolation
Xq = xbins(1):binsize_final:xbins(end); % after interpolation


% ================ SAVE ORIGINAL ERSION
OUTSTRUCT_XCOV.XcovgramBase_orig = OUTSTRUCT_XCOV.XcovgramBase;
OUTSTRUCT_XCOV.XcovgramWN_orig = OUTSTRUCT_XCOV.XcovgramWN;
%%

for i=1:length(OUTSTRUCT_XCOV.bnum)
    xgram_base = OUTSTRUCT_XCOV.XcovgramBase{i};
    xgram_wn = OUTSTRUCT_XCOV.XcovgramWN{i};
    
    
    % ==================== 1) interpolate xcov traces
    % INTERP PARAMS.
    
    % 1) INTERPOLATE BEFORE SMOOTHING
    xgram_base = interp1(xbins, xgram_base', Xq)';
    xgram_wn = interp1(xbins, xgram_wn', Xq)';
    
    % =================== 2) SMOOTH
    % SMOOTH PARAMS (gaussian)
    N = 2*dosmooth_sigma/binsize_final;
    win = gausswin(N);
    win = win./sum(win);
    
    % DO SMOOTH
    assert(size(xgram_wn,1)==size(xgram_base,1));
    for kkk=1:size(xgram_base,1)
        xgram_base(kkk,:) = conv(xgram_base(kkk, :), win, 'same');
        xgram_wn(kkk,:) = conv(xgram_wn(kkk, :), win, 'same');
    end
    
    % ====== DOWNSAMPLE IF BINSIZE IS <1ms
    if binsize_final<0.001
        Xq_final = Xq(1):0.001:Xq(end);
        xgram_base = interp1(Xq, xgram_base', Xq_final)';
        xgram_wn = interp1(Xq, xgram_wn', Xq_final)';
    end
    
    % ================== OUTPUT
    OUTSTRUCT_XCOV.XcovgramBase{i} = xgram_base;
    OUTSTRUCT_XCOV.XcovgramWN{i} = xgram_wn;
end
%%

PARAMS.Xcov_ccLags_orig = PARAMS.Xcov_ccLags;
PARAMS.Xcov_ccLags = Xq_final;