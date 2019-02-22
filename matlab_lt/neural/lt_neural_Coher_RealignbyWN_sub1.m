function [ytmp, xtmp] = lt_neural_Coher_RealignbyWN_sub1(OUTSTRUCT, PARAMS, fieldtoget)
%% outputs matrices that you should break up into cells if desired.
tbins = PARAMS.tbins;
binsize = tbins(2)-tbins(1);

% t1 = -round(min(OUTSTRUCT.WNonset)/binsize);
t2 = round(max(OUTSTRUCT.WNonset)/binsize); % max number of shift

% ytmp = nan(length(PARAMS.tbins)+t1+t2, length(PARAMS.ffbins), length(OUTSTRUCT.bnum));
ytmp = nan(length(PARAMS.tbins)+t2, length(PARAMS.ffbins), length(OUTSTRUCT.bnum));
xdefault = [1:length(PARAMS.tbins)]+t2;

%% goes thru all cases

for i=1:length(OUTSTRUCT.bnum)
    
%    OUTSTRUCT.PhiMat;
%    OUTSTRUCT.CohMean_WNminusBase;
%    OUTSTRUCT.CohMean_Base;
%    OUTSTRUCT.CohMean_WN;

   ymat = OUTSTRUCT.(fieldtoget){i};
%    ymat_wn = OUTSTRUCT.(fieldtoget){i};
   
   wntime = OUTSTRUCT.WNonset(i);
   xthis = xdefault - round(wntime/binsize); % shift time base depending on time of WN.
   
   if ~isempty(ymat)
   ytmp(xthis, :, i) = ymat;
   end
   
%    tbinsnew = PARAMS.tbins - wntime;
%    
%    tshift = round(wntime/binsize);
   
   
   
   % -- round to nearest bin
   
   
   %    OUTSTRUCT.Spec1Mean_Base
%    OUTSTRUCT.S
    
end

%% 
% ----- CONFIRM THAT ALL TIME BINS FILLED. not necessary I guess...
if (0)
    tmp = squeeze(ytmp(:,1,:));
assert(~any(all(isnan(tmp')))); % i.e. all time bins are filled
end

% --- only keep time bins that have data for all cases
indsrows = ~cellfun(@(x)isempty(x), OUTSTRUCT.(fieldtoget)); % these are the case that actually had data...
tmp = squeeze(ytmp(:,1,:));
indstoremove = any(isnan(tmp(:,indsrows)')); % any time bins that have at least one nan...
ytmp(indstoremove, :,:) = [];

%% ==== figure out new time bins.
% xtmp = 1:max(xdefault);
% xtmp(xdefault) = PARAMS.tbins;
% xtmp(indstoremove) = [];

xtmp = nan(1, max(xdefault));
assert(length(xtmp)==length(indstoremove));
xtmp(xdefault) = PARAMS.tbins;
xtmp(indstoremove) = [];

assert(~isnan(xtmp(end)));
x1 = xtmp(end)-binsize*(length(xtmp)-1);
xend = xtmp(end);
xtmp = single(linspace(x1, xend, length(xtmp)))

% === convert from mat to cell
ytmpcell = cell(size(ytmp,3),1);
for j=1:size(ytmp,3)
   ytmpcell{j} = ytmp(:,:,j);
end
ytmp = ytmpcell;
