function lt_neural_POPLEARN_XCov_ReAlign(OUTSTRUCT_XCOV, OUTSTRUCT, SwitchCohStruct, ...
    SwitchStruct, PARAMS, useallwn)
%% lt 2/15/19 -


%% outputs matrices that you should break up into cells if desired.
tbins = PARAMS.xcenters_gram;
binsize = tbins(2)-tbins(1);

t1 = -round(min(OUTSTRUCT.WNonset)/binsize);
t2 = round(max(OUTSTRUCT.WNonset)/binsize); % max number of shift

ytmp = nan(length(tbins)+t1+t2, length(PARAMS.Xcov_ccLags), length(OUTSTRUCT_XCOV.bnum));
xdefault = [1:length(tbins)]+t2;


for i=1:length(OUTSTRUCT_XCOV.bnum)
    
    % ========= find timing.
    bnum = OUTSTRUCT_XCOV.bnum(i);
    enum = OUTSTRUCT_XCOV.enum(i);
    sw = OUTSTRUCT_XCOV.switch(i);
    
    indstmp = OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & OUTSTRUCT.switch==sw;
    wntime = unique(OUTSTRUCT.WNonset(indstmp));
    assert(length(wntime)==1);
    
    
    % =========
    xmat_base = OUTSTRUCT_XCOV.XcovgramBase{i};
    xmat_wn = OUTSTRUCT_XCOV.XcovgramWN{i};
    
    
    % ==========
    xthis = xdefault - round(wntime/binsize); % shift time base depending on time of WN.
    
    ytmp(xthis, :, i) = xmat_base;
    
    
    %    tbinsnew = PARAMS.tbins - wntime;
    %
    %    tshift = round(wntime/binsize);
    
    
    
    % -- round to nearest bin
    
    
    %    OUTSTRUCT.Spec1Mean_Base
    %    OUTSTRUCT.S
    
    
end

%% =========== continue modifying code from here. stopped becuase did not extract enough close to WN.
% ----- CONFIRM THAT ALL TIME BINS FILLED. not necessary I guess...
tmp = squeeze(ytmp(:,1,:));
assert(~any(all(isnan(tmp')))); % i.e. all time bins are filled

% --- only keep time bins that have data for all cases
tmp = squeeze(ytmp(:,1,:));
indstoremove = any(isnan(tmp')); % any time bins that have at least one nan...
ytmp(indstoremove, :,:) = [];

% ==== figure out new time bins.
xtmp = 1:max(xdefault);
xtmp(xdefault) = PARAMS.tbins;
xtmp(indstoremove) = [];

% === convert from mat to cell
ytmpcell = cell(size(ytmp,3),1);
for j=1:size(ytmp,3)
    ytmpcell{j} = ytmp(:,:,j);
end
ytmp = ytmpcell;