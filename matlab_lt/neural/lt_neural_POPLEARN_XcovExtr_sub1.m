function [xcovgram_base, xcovgram_wn] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
    dosmooth, dosmooth_sigma, inds_base, inds_WN, PARAMS, xcovver, datbase, ...
    np, hilosplit_shuffver)

if ~exist('hilosplit_shuffver', 'var')
    hilosplit_shuffver = 1;
end
%% lt - outputs xcovgram

nwinds = size(datthis.ccRealAllPair_allwind,2);
xcovgram_base = nan(nwinds, length(datbase)); % win x lags
xcovgram_wn= nan(nwinds, length(datbase));
for ww=1:nwinds
    
    datmat_real = datthis.ccRealAllPair_allwind{np, ww};
    datmat_shuff = datthis.ccShiftAllPair_allwind{np, ww};
    
    if strcmp(xcovver, 'coherency')
        % ---- autocovariance
        datmat_auto1 = datthis.ccRealAllPair_Auto1_allwind{np, ww};
        datmat_auto2 = datthis.ccRealAllPair_Auto2_allwind{np, ww};
        datcell_auto_real = {datmat_auto1, datmat_auto2};
        
        % --- autocovariance (SHIFTED)
        datmat_auto1 = datthis.ccShiftAllPair_Auto1_allwind{np, ww};
        datmat_auto2 = datthis.ccShiftAllPair_Auto2_allwind{np, ww};
        datcell_auto_shift = {datmat_auto1, datmat_auto2};
    else
        datcell_auto_real = {};
        datcell_auto_shift = {};
    end
    [db, dw] = ...
        lt_neural_POPLEARN_XCov_sub1(datmat_real, datmat_shuff, dosmooth, ...
        dosmooth_sigma, inds_base, inds_WN, PARAMS.Xcov_ccLags, 0, ...
        xcovver, datcell_auto_real, datcell_auto_shift, hilosplit_shuffver);
    
    xcovgram_base(ww, :) = db;
    xcovgram_wn(ww, :) = dw;
    
    if any(imag(db)~=0)
        keyboard
    end
end

% ======= make single
xcovgram_base = single(xcovgram_base);
xcovgram_wn = single(xcovgram_wn);


