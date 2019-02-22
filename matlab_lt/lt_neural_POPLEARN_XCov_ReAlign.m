function [OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XCov_ReAlign(OUTSTRUCT_XCOV, OUTSTRUCT, PARAMS)
%% lt 2/15/19 -


%% outputs matrices that you should break up into cells if desired.
tbins = PARAMS.xcenters_gram;
binsize = tbins(2)-tbins(1);

% t1 = -round(min(OUTSTRUCT.WNonset)/binsize);
t2 = round(max(OUTSTRUCT.WNonset)/binsize); % max number of shift

% ytmp = nan(length(tbins)+t1+t2, length(PARAMS.Xcov_ccLags), length(OUTSTRUCT_XCOV.bnum));
% ytmp_wn = nan(length(tbins)+t1+t2, length(PARAMS.Xcov_ccLags), length(OUTSTRUCT_XCOV.bnum));
ytmp = nan(length(tbins)+t2, length(PARAMS.Xcov_ccLags), length(OUTSTRUCT_XCOV.bnum));
ytmp_wn = nan(length(tbins)+t2, length(PARAMS.Xcov_ccLags), length(OUTSTRUCT_XCOV.bnum));
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
%     if wntime>0.05
%         keyboard
%     end
    xthis = xdefault - round(wntime/binsize); % shift time base depending on time of WN.
    
    ytmp(xthis, :, i) = xmat_base;
    ytmp_wn(xthis, :, i) = xmat_wn;
    
    
    %    tbinsnew = PARAMS.tbins - wntime;
    %
    %    tshift = round(wntime/binsize);
    
    
    
    % -- round to nearest bin
    
    
    %    OUTSTRUCT.Spec1Mean_Base
    %    OUTSTRUCT.S
    
    
end

%% =========== continue modifying code from here. stopped becuase did not extract enough close to WN.
% ----- CONFIRM THAT ALL TIME BINS FILLED. not necessary I guess...
if (0)
tmp = squeeze(ytmp(:,1,:));
assert(~any(all(isnan(tmp')))); % i.e. all time bins are filled
end

% --- only keep time bins that have data for all cases
tmp = squeeze(ytmp(:,1,:));
indstoremove = any(isnan(tmp')); % any time bins that have at least one nan...
% -- BASE AND WN
ytmp(indstoremove, :,:) = [];
ytmp_wn(indstoremove, :,:) = [];

% === convert from mat to cell
ytmpcell = cell(size(ytmp,3),1);
for j=1:size(ytmp,3)
    ytmpcell{j} = ytmp(:,:,j);
end
ytmp = ytmpcell;

% === convert from mat to cell
ytmpcell = cell(size(ytmp_wn,3),1);
for j=1:size(ytmp_wn,3)
    ytmpcell{j} = ytmp_wn(:,:,j);
end
ytmp_wn = ytmpcell;

% ===================== add to OUTSTDUCT
OUTSTRUCT_XCOV.XcovgramBase_alignWN = ytmp;
OUTSTRUCT_XCOV.XcovgramWN_alignWN = ytmp_wn;

%% ==== figure out new time bins.
% --- assumes that no WN starts after syl onset
assert(all(OUTSTRUCT.WNonset>0));

xtmp = nan(1, max(xdefault));
assert(length(xtmp)==length(indstoremove));
xtmp(xdefault) = PARAMS.xcenters_gram;
xtmp(indstoremove) = [];

% --- fill in the nan,...
assert(~isnan(xtmp(end)));
x1 = xtmp(end)-binsize*(length(xtmp)-1);
xend = xtmp(end);
xthis = single(linspace(x1, xend, length(xtmp)));
PARAMS.xcenters_gram_alignWN = xthis;
