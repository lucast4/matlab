function [segextract, xtimes] = lt_neural_QUICK_SpkBinned(segextract, ...
    maxdur, binsize, convertosingle, dojitter, jitterwindSize)
%% lt 1/17/18 - takes spktime and outputs binned (1ms res)

if ~exist('dojitter', 'var')
    dojitter = 0;
end

if ~exist('binsize', 'var')
    binsize = 0.001;
end

if ~exist('maxdur', 'var')
    maxdur = []; % in sec, max dur to get bins (i.e. for xbins)
end

%% NOTE if TIME WARPED:
% 1) bins start at same time across trials (i.e. alignment point - motif predur). throws out spikes that occur before this
% 2) bins end at last spike across trials - not perfect, since there might be silence past that ...

%% -- run

if isempty(maxdur)
    % --- what is maximum common trial dur?
    stimes = {segextract.spk_Times};
    maxtime = cellfun(@max, stimes);
    maxtime = max(maxtime);
else
    maxtime = maxdur;
end

% --- initiate xbins
% xedges = 0:0.001:ceil(maxtime*1000)/1000;
% maxtime-mod(maxtime, binsize)
xedges = 0:binsize:maxtime;
xtimes = (xedges(2)-binsize/2):binsize:(xedges(end)-binsize/2);

%% do jitter?
if dojitter==1
    spktimes = {segextract.spk_Times};
   
   % -- convert to one long array
   tnumall = [];
   stimeall = [];
   clustnumall = [];
   for j=1:length(spktimes)
           stimeall = [stimeall spktimes{j}];
           tnumall = [tnumall j*ones(size(spktimes{j}))];
           clustnumall = [clustnumall segextract(j).spk_Clust];
   end
   
   % ======= for each jitter window perform jitter
   jitterWindOff = jitterwindSize:jitterwindSize:maxdur;
   if maxdur-jitterWindOff(end)>jitterwindSize/2
       % -- take the last window
       jitterWindOff = [jitterWindOff jitterWindOff(end)+jitterwindSize];
   elseif maxdur-jitterWindOff(end)==0
       % -- do nothing       
   else
       % -- combine with last window
       jitterWindOff(end) = maxdur;
   end
   jitterWindOn = jitterWindOff-jitterwindSize;
   
   % ============ got thru each jitter window.
   for j=1:length(jitterWindOn)
      jitOn = jitterWindOn(j);
      jitOff = jitterWindOff(j);
      
      % ============ COLLECT ALL SPIKES IN THIS WINDOW
      sInds = find(stimeall>jitOn & stimeall<=jitOff);
      
      sTimes = stimeall(sInds);
      cNums = clustnumall(sInds);
      
      % =========== mix up the spike times and put them back.
      indrand = randperm(length(sTimes));
      stimeall(sInds) = sTimes(indrand);
      clustnumall(sInds) = cNums(indrand);
   end
   
   % ======== PUT SPIKES BACK INTO ORIGINAL SEGEXTRACT
   for j=1:length(segextract)
       
       Sthis = stimeall(tnumall == j);
       Cthis = clustnumall(tnumall ==j);
       
       [~, indsort] = sort(Sthis);
       
       Sthis = Sthis(indsort);
       Cthis = Cthis(indsort);
       
       segextract(j).spk_Times = Sthis;
       segextract(j).spk_Clust = Cthis;
%        assert(length(unique(segextract(j).spk_Clust))==1, 'then the mapping between clust and spike will be messed up..');
   end
end

%% for each trial get counts
ntrials = length(segextract);
SpkCounts = nan(ntrials, length(xedges)-1);
for t =1:ntrials
    y = histc(segextract(t).spk_Times, xedges);
    y = y(1:end-1);
    
    if convertosingle==1
        y = single(y);
    else
        y = int8(y);
    end
    
    SpkCounts(t,:) = y;
    segextract(t).spk_Binned = y';
    segextract(t).spk_Binned_x = xtimes';
end

tmp = isnan(SpkCounts);
assert(~any(isnan(tmp(:))), 'asdfas');

%%

