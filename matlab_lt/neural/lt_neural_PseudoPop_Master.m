close all; clear MOTIFSTATS_Compiled;
lt_neural_ExtractMotifs_Regular;

%% ========== FOR EACH SYLLABLE, GET POPULATION VECTOR (VECTOR OF MEANS)
close all

NormToNeurMean = 1;
% premotorWind = [-0.06 -0.02];
% premotorWind = [-0.03 0.01];
premotorWind = [-0.03 0.02];
plotnans = 1;
birdnum = 1;

[FRmatMotifByNeur, AllMotifRegexp, AllNeurLocation] = ...
    lt_neural_v2_PseudoPop_Master(MOTIFSTATS_Compiled, NormToNeurMean, ...
    premotorWind, plotnans, birdnum);

%% ================ SAME-TYPE VS. DIFF SUMMARY

plotOn=1;
thisloc = 'RA';
distmetric = 'euclidean';
distmetric = 'correlation';

[dp, FRmatDist] = lt_neural_v2_PseudoPop_SameDiff(FRmatMotifByNeur, AllMotifRegexp, ...
    AllNeurLocation, thisloc, plotOn, distmetric);



%% ############################# SYL ENCODING AS FUNCTION OF WINDOW
if (1)
    WindList = [-0.1:0.01:0.08; -0.06:0.01:0.12]';
    DprimeAll = [];
    thisloc = 'RA';
    distmetric = 'euclidean';
    distmetric = 'correlation';
    
    for j=1:size(WindList,1)
        windthis = WindList(j,:);
        
        % ==================== get pop vector
        [FRmatMotifByNeur, AllMotifRegexp, AllNeurLocation] = ...
            lt_neural_v2_PseudoPop_Master(MOTIFSTATS_Compiled, 1, windthis, 0);
        
        % ====================== get dprime (diff minus same)
        plotOn=0;
        
        dp = lt_neural_v2_PseudoPop_SameDiff(FRmatMotifByNeur, AllMotifRegexp, ...
            AllNeurLocation, thisloc, plotOn, distmetric);
        
        DprimeAll = [DprimeAll; dp];
    end
    
    % ============== PLOT
    lt_figure; hold on;
    title(thisloc);
    xlabel('wind center');
    ylabel('dprime(diff - same');
    
    plot(mean(WindList,2), DprimeAll, '-ok');
    
end
%% ================ VARIOUS RAW PLOTS

lt_neural_v2_PseudoPop_PLOT;


%% ================ HEIRARCHICAL CLUSTERING
close all;
distmetric = 'euclidean';
distmetric = 'correlation';
locthis = 'RA';
inds = find(strcmp(AllNeurLocation, locthis));

Z = linkage(FRmatMotifByNeur(:,inds), 'single', distmetric);
lt_figure; hold on;
dendrogram(Z, 'Labels', AllMotifRegexp);
% dendrogram(Z);
rotateXLabels(gca, 90);

tmp = get(gca,'XTickLabel');

% ---- calcualte cophenetic correlation (i..e how well does clustering
% represent original distances?)

Y = pdist(FRmatMotifByNeur, distmetric);
c = cophenet(Z, Y);

%% =============== PICK ONE SYL, PLOT ITS DISTANCE TO EVERY OTHER SYL
seedmotif = 'aa(b)hh';
locthis = 'RA';
distmetric = 'correlation';

indneur = strcmp(AllNeurLocation, locthis);
indmotif = strcmp(AllMotifRegexp, seedmotif);

frmat = FRmatMotifByNeur(:,indneur);

if strcmp(distmetric, 'correlation')
   
    % === get correlation between all motifs
    rhomat = corr(frmat');
    corrvalues = rhomat(indmotif, :);
    
    lt_figure; hold on;
    lt_plot_bar(1:length(corrvalues), corrvalues);
    
    set(gca, 'XTick', 1:length(AllMotifRegexp), 'XTickLabel', AllMotifRegexp);
    title([locthis ', pop corr to ' seedmotif]);
end


%% ===== CORRELATION BETWEEN DISTANCE AND GENERALIZATION LEARNING
% ======================== 1) SAVE DISTANCE MATRIX FOR EACH BIRD


% ======================== 2)

%% #########################################################
%% ############################## RUNNING POPULATION ANALYSES

birdnum = 1;
norm_zscore = 1; % -- 0: no norm (just nspks) 1: zscore each neuron (across all time bins)
    % if 2, then norms by min and max FR (i.e. 0 to 1); 
    % 3: zscore separately for each neuron and time bin (across motifs) --

binsize = 0.01; % for spikes
WindToTake = [-0.1 0.05]; % relative to syl onset

wh44tempfix = 0;
if wh44tempfix==1
    % -- then add some motifs, these are the 3 motifs that are convergent
    % while controlling for following syllable
    motifstoadd = {'j(j)j', 'b(j)j', 'h(m)d', 'm(m)d', 'b(n)h', 'j(n)h'};
end

RemoveNeuronsWithoutAllMotifs=1; % [default = 0] if 0 then first removes neurons
% that lack lot of motifs. then removes motifs that lack lots of neurons.


%% =====
% =======================================


close all;
assert(length(SummaryStruct.birds) == length(MOTIFSTATS_Compiled.birds), 'asfsd');

i=birdnum;

assert(length(SummaryStruct.birds(i).neurons) == length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons), 'asfd');

% ====== for this bird
Motifstats = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS;

motiflist = Motifstats.params.motif_regexpr_str;

% ===== temporary fix
birdname = SummaryStruct.birds(i).birdname;
if strcmp(birdname, 'wh44wh39') & wh44tempfix==1
   % === append motifs
   motiflist = [motiflist motifstoadd];
end
    
    
% -- get list of single syls
[~, ~, ~, SingleSyls] = ...
    lt_neural_v2_extractSameType(motiflist, {motiflist{1}});

singlesylUnique = unique(SingleSyls);
numneurons = length(Motifstats.neurons);

% ======== FIRST, COLLECT BRAIN REGIONS FOR NEURONS
NeurBrainRegions = {};
for nn=1:numneurons
    NeurBrainRegions = [NeurBrainRegions SummaryStruct.birds(i).neurons(nn).NOTE_Location];
end



%% ======================================== OLD VERSION, ONLY GET SAME TYPES

if (1)

ClassStruct = struct;
cc= 0 ;
NNueronsCount = [];
for j=1:length(singlesylUnique)
    
    ssyl = singlesylUnique{j};
    
    % ---- figure out how many contexts for this
    indsmotifs = find(strcmp(SingleSyls, ssyl));
    
    if length(indsmotifs)==1
        % then only one context ...
        continue
    end
    
    % ======= for each motif collect population vector at each time point
    % ---- bin spike params
    maxdur = Motifstats.params.motif_predur+Motifstats.params.motif_postdur;
%     maxdur = WindToTake(2) - WindToTake(1);
    converttosingle = 1;
    numtbins = floor((WindToTake(2) - WindToTake(1))/binsize);
    
    NspkNeurTimeMotif = nan(numneurons, numtbins, length(indsmotifs)); % neurons x timebins x motif
    for k = 1:length(indsmotifs)
        
        mm = indsmotifs(k);
        % ================== FOR THIS MOTIF, GO THRU ALL NEURONS, FOR EACH
        % NEURONS COLLECT A FR VECTOR OVER TIME. COMBINE ACROSS ALL NEURONS
        % TO FINALIZE
        
        
        for nn=1:numneurons
            
            segextract = Motifstats.neurons(nn).motif(mm).SegmentsExtract;
            if isempty(segextract)
                continue
            end
               
            [segextract, xtimes] = lt_neural_QUICK_SpkBinned(segextract, ...
                maxdur, binsize, converttosingle);
            
            Nspkall = [segextract.spk_Binned];
            NspkMeanVect = mean(Nspkall,2);
            
            % --- take desired FR window
            indtmp = xtimes >= Motifstats.params.motif_predur+WindToTake(1) ...
                & xtimes <= Motifstats.params.motif_predur+WindToTake(2);
            NspkMeanVect = NspkMeanVect(indtmp);
            
            NspkNeurTimeMotif(nn, :, k) = NspkMeanVect;
            %             NspkNeurTimeMat = [NspkNeurTimeMat; NspkMeanVect'];
        end
        
        
    end
%     assert(~any(isnan(NspkNeurTimeMotif(:))), 'sdfasd');
    
    % =================== OUTPUT
    cc = cc+1;
    ClassStruct.branch(cc).singlesyl = ssyl;
    ClassStruct.branch(cc).indsmotifs = indsmotifs;
    ClassStruct.branch(cc).motifnames = motiflist(indsmotifs);
    ClassStruct.branch(cc).NspkNeurTimeMotif = NspkNeurTimeMotif;
    
   
end
ClassStruct.NeurBrainRegions = NeurBrainRegions;
end

%% ============================ COLLECT FR vs TIME mat for all motifs
% ClassStruct = struct;
% cc= 0 ;
% NNueronsCount = [];
i=birdnum;

Motifstats = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS;
% ---- bin spike params
maxdur = Motifstats.params.motif_predur+Motifstats.params.motif_postdur;
converttosingle = 1;
    numtbins = floor((WindToTake(2) - WindToTake(1))/binsize);

NspkNeurTimeMotif = nan(numneurons, numtbins, length(motiflist)); % neurons x timebins x motif
for mm=1:length(motiflist)
    
    % ======= for each motif collect population vector at each time point
    
    for nn=1:numneurons
        
        segextract = Motifstats.neurons(nn).motif(mm).SegmentsExtract;
        
        if isempty(segextract)
            continue
        end
        
        [segextract, xtimes] = lt_neural_QUICK_SpkBinned(segextract, ...
            maxdur, binsize, converttosingle);
        
        Nspkall = [segextract.spk_Binned];
        NspkMeanVect = mean(Nspkall,2);
        
        % -- take desired FR window
        indtmp = xtimes >= Motifstats.params.motif_predur+WindToTake(1) ...
                & xtimes <= Motifstats.params.motif_predur+WindToTake(2);
         NspkMeanVect = NspkMeanVect(indtmp);

        
        NspkNeurTimeMotif(nn, :, mm) = NspkMeanVect;
        %             NspkNeurTimeMat = [NspkNeurTimeMat; NspkMeanVect'];
    end
    
end
% ClassStruct.NeurBrainRegions = NeurBrainRegions;

% assert(~any(isnan(NspkNeurTimeMotif(:))), 'sdfasd');

% ================= THROW OUT NEURONS THAT LACK A LOT OF DATA
neuronstoremove = [];
for nn=1:size(NspkNeurTimeMotif,1)
    
    tmp = (NspkNeurTimeMotif(nn,:,:));
    disp(['nans: ' num2str(sum(isnan(tmp(:)))) '/' num2str(length(tmp(:)))]);
    
    % --- if a neuron is missing more than 3/4 of data, then remove
    numbermissing =  sum(isnan(tmp(:)));
    numbertot = length(tmp(:));
    
    if RemoveNeuronsWithoutAllMotifs==1
        threshval = 0;
    else
        threshval = 0.75;
    end
    if numbermissing/numbertot > threshval
        neuronstoremove = [neuronstoremove nn];
    end
    
end
% -------- remove those neurons
NeuronsRemaining = 1:size(NspkNeurTimeMotif,1);
if ~isempty(neuronstoremove)
disp(['[lack >75% data] REMOVING neurons ' num2str(neuronstoremove)]);
pause
NspkNeurTimeMotif(neuronstoremove, :,:) = [];
NeurBrainRegions(neuronstoremove) = [];
NeuronsRemaining(neuronstoremove) = [];
end


% ========================================= THROW OUT MOTIFS THAT LACK DATA
motiftoremove = [];
for mm=1:length(motiflist)
% lt_figure; hold on;
% title(num2str(mm));
% spy(NspkNeurTimeMotif(:,:,mm));

if any(isnan(NspkNeurTimeMotif(:,:,mm)))
    motiftoremove = [motiftoremove mm];
end

end

if ~isempty(motiftoremove)
   
    disp(['NOTE: removing motifs ' num2str(motiftoremove) ' since lack data']);
     pause;
     
     MotiflistAfterRemove = motiflist;
     MotiflistAfterRemove(motiftoremove) = [];
     NspkNeurTimeMotif(:,:,motiftoremove) = [];
else
     MotiflistAfterRemove = motiflist;
end

%% reset params

numneurons = size(NspkNeurTimeMotif,1);
numtbins = size(NspkNeurTimeMotif,2);

%% ============ perform normalization
if norm_zscore==1
   % -- zscore norm
    for nn = 1:numneurons
        
        % -- get mean and std
        tmp = NspkNeurTimeMotif(nn, :,:);
        ymean = mean(tmp(:));
        ystd = std(tmp(:));
        
        % --- normlaize
        tmp = squeeze(tmp);
        NspkNeurTimeMotif(nn, :,:) = (tmp - ymean)./ystd; % put back into mat
        
    end
elseif norm_zscore==2
    % --- min-max normalization
    for nn = 1:numneurons
        
        % -- get mean and std
        tmp = NspkNeurTimeMotif(nn, :,:);
        ymin = prctile(tmp(:), 2.5); % 2.5 percentile
        ymax = prctile(tmp(:), 97.5);

        % --- normlaize
        tmp = squeeze(tmp);
        tmp = (tmp - ymin)./(ymax-ymin);
        
        NspkNeurTimeMotif(nn, :,:) = tmp;  % put back into mat
        
    end
elseif norm_zscore==3
    % --- for each neuron and time bin, normalize based on firing across
    % motifs
    for nn=1:numneurons
       for tt = 1:numtbins
        
           % --- get mean and std
           tmp = NspkNeurTimeMotif(nn, tt, :);
           ymean = mean(tmp(:));
           ystd = std(tmp(:));
           
           % -- nromalize
           tmp = squeeze(tmp);
           NspkNeurTimeMotif(nn, tt, :) = (tmp - ymean)./ystd;
           
       end
    end
end

%% ==== perform normalization
if (1) % -- this is old version ...
if norm_zscore==1
% --- each neuron, normalize using all time bins across all motifs
    numbranch = length(ClassStruct.branch);
    
for nn=1:numneurons
   
    for b = 1:numbranch
        
       tmp = squeeze(ClassStruct.branch(b).NspkNeurTimeMotif(nn, :, :));
       tmp = tmp(:);
       
       ymean = mean(tmp);
       ystd = std(tmp);
       
      
       nummotifs = length(ClassStruct.branch(b).motifnames);
       for mm=1:nummotifs
          
           tmp = ClassStruct.branch(b).NspkNeurTimeMotif(nn, :, mm);
           tmp = (tmp-ymean)./ystd;
           
 ClassStruct.branch(b).NspkNeurTimeMotif(nn, :, mm) = tmp;
           
       end
       
    end
   
end


% --- I know all branches extracted must have the same number of neurons
% currently
end
end



%% ###################### PLOT
% close all;
numbranches = length(ClassStruct.branch);
oldversion = 0; % only diff -  normlaizes only using motifs that are in branches. [LEAVE at 0]
usecorr = 1; % if 0, then uses euclidian. 1 = corr; 2 = dot product, 3: cosine angle
dorobustcorr =1; % if 1, then 

% ======================= 1) PLOT POPULATION OVER TIME (HEAT MAP
bregion = 'RA';
indsneur = strcmp(NeurBrainRegions, bregion);

for j=1:numbranches
    lt_figure; hold on;
    
    numclasses = length(ClassStruct.branch(j).motifnames);
    for jj=1:numclasses
        motif = ClassStruct.branch(j).motifnames{jj};
        lt_subplot(4,2,jj); hold on;
        title([motif '-' bregion]);
        xlabel('tbin');
        ylabel('neuron');
        
        if oldversion==1
            nspkmat = ClassStruct.branch(j).NspkNeurTimeMotif(indsneur,:, jj);
        else
            indtmp = strcmp(MotiflistAfterRemove, motif);
            nspkmat = NspkNeurTimeMotif(indsneur,:, indtmp);
        end
        n = size(nspkmat, 1);
        t = size(nspkmat, 2);
        
        imagesc(1:t, 1:n, nspkmat);
        axis tight;
        colorbar
    end
    
    % ============ mean FR across neurons
    lt_subplot(4,2, numclasses+1); hold on ;
    ylabel(['mean FR across neurons']);
    
    for k=1:numclasses
        if oldversion==1
            frmean = mean(ClassStruct.branch(j).NspkNeurTimeMotif(indsneur, :, k),1);
            
        else
            motif = ClassStruct.branch(j).motifnames{k};
            indtmp = strcmp(MotiflistAfterRemove, motif);
            frmean = mean(NspkNeurTimeMotif(indsneur,:, indtmp),1);
            
        end
        plot(frmean, '-o', 'Color', [rand rand rand]);
    end
    
    
    % ================ corelation between population vectors
    lt_subplot(4,2,numclasses+2); hold on;
    ylabel('pairwise pop correlation between motifs');
    for k=1:numclasses
        for kk=k+1:numclasses
            
            if oldversion==1
                frmat1 = ClassStruct.branch(j).NspkNeurTimeMotif(indsneur, :, k);
                frmat2 = ClassStruct.branch(j).NspkNeurTimeMotif(indsneur, :, kk);
            else
                motif1 = ClassStruct.branch(j).motifnames{k};
                motif2 = ClassStruct.branch(j).motifnames{kk};
                
                indtmp1 = strcmp(MotiflistAfterRemove, motif1);
                indtmp2 = strcmp(MotiflistAfterRemove, motif2);
                
                frmat1 = NspkNeurTimeMotif(indsneur, :, indtmp1);
                frmat2 = NspkNeurTimeMotif(indsneur, :, indtmp2);
                
                neuronIDlist = NeuronsRemaining(indsneur);
            end
            
            if usecorr==1
                % =============== 1) use correlation metric
                rho = corr(frmat1, frmat2);
                y = diag(rho);
                ylim([-1 1]);
                
%                 % ================ robust correlation (leave one out)
%                 y1 = frmat1(:, 5);
%                 y2 = frmat2(:, 5);
%                 figure; hold on;
%                 plot(y1, '-ok');
%                 plot(y2, '-ob');
            elseif usecorr == 0
                % ====== 2) use euclidian distance
                tmp = frmat1 - frmat2;
                y = sqrt(sum(tmp.^2, 1));
            elseif usecorr ==2
                % ===== 3) use dot product.
                y = sum(frmat1.*frmat2,1);
            elseif usecorr ==3
                % ====4) norm dot product (i.e. cosine angle)
                y = sum(frmat1.*frmat2,1);
                norm1 = sqrt(sum(frmat1.^2,1));
                norm2 = sqrt(sum(frmat2.^2,1));
                y = y./(norm1.*norm2);
            end
            plot(1:length(y), y, '-ok');
        end
    end
    lt_plot_zeroline;
    
    % ================= GET AUTOCORRELATION OF POPULATION ACTIVITY
    for k=1:numclasses
            motif = ClassStruct.branch(j).motifnames{k};
            indtmp = strcmp(MotiflistAfterRemove, motif);
            frmat = NspkNeurTimeMotif(indsneur,:, indtmp);
            
            rhomat = corr(frmat);
            lt_subplot(4,2,numclasses+3+k-1); hold on;
            imagesc(1:t, 1:t, rhomat);
   axis tight
        
    end
    
    
end

% ================= PLOT DIFFERENT TIME SLICES
timesliceinterval = 3;
for j=1:numbranches
    lt_figure; hold on;
    
    numclasses = length(ClassStruct.branch(j).motifnames);
    pcols = lt_make_plot_colors(numclasses, 0, 0);
    % ====== go thru each time slice
    timeslices = 1:timesliceinterval:size(NspkNeurTimeMotif,2);
    for tt = 1:length(timeslices)
       timebin = timeslices(tt);
        lt_subplot(3,4,tt); hold on;
        title([bregion ',timebin' num2str(timebin)]);
       
       % ====== plot all classes for this times lice
    for jj=1:numclasses
        motif = ClassStruct.branch(j).motifnames{jj};
%         title([motif '-' bregion]);
%         xlabel('tbin');
%         ylabel('neuron');
        
            indtmp = strcmp(MotiflistAfterRemove, motif);
            nspkmat = NspkNeurTimeMotif(indsneur,timebin, indtmp);

        plot(1:length(nspkmat), nspkmat, '-o', 'Color', pcols{jj});
        
    end   
        
        legend(gca, ClassStruct.branch(j).motifnames)
    end
end
%% ==================== RUNNING CORREALTION BETWEEN ALL DIFF PAIRS

BrainRegionList = {'LMAN', 'RA'};

RhoStruct = struct;
count =1; % to keep track of data points (e.g. a given motif pair can have multiple brain regions)
motifpairID = 0; % to keep track of motif pairs
% ==== go thru all pairs, collect running corr, and not if is same or diff
for j=1:length(MotiflistAfterRemove)
    for jj=j+1:length(MotiflistAfterRemove)
        
        % ============== 1) same or diff
        motif1 = MotiflistAfterRemove{j};
        motif2 = MotiflistAfterRemove{jj};
        issame = lt_neural_QUICK_SylPairType(motif1, motif2);
        motifpairID = motifpairID+1;

        % -- what is single syl?
        
        
        % ============= 2) collect correlation [separately by brain region]
        for bb = 1:length(BrainRegionList)
            bregionthis = BrainRegionList{bb};
            
            indsneur = strcmp(NeurBrainRegions, bregionthis);
            
            indmotif1 = strcmp(MotiflistAfterRemove, motif1);
            indmotif2 = strcmp(MotiflistAfterRemove, motif2);
            
            frmat1 = NspkNeurTimeMotif(indsneur, :, indmotif1);
            frmat2 = NspkNeurTimeMotif(indsneur, :, indmotif2);
            
            % -------- calculate rho
            rho = corr(frmat1, frmat2);
            y = diag(rho);
            
            RhoStruct.motifpair(count).rhovstime = y;
            RhoStruct.motifpair(count).motif1 = motif1;
            RhoStruct.motifpair(count).motif2 = motif2;
            RhoStruct.motifpair(count).issame = issame;
            RhoStruct.motifpair(count).motifpairID = motifpairID;
            RhoStruct.motifpair(count).brainregion = bregionthis;
            
            count = count+1;
        end
    end
end

%% ========== ASK HOW WELL ACTIVITY CORRELATES, FOR EACH TIME POINT, ABOVE RANDOM PAIRS

% --- for each motif pair, get all cases of [motif1 vs. rest] and [motif2
% vs. rest] and take average of that. substract that from motif1 vs.
% motif2.

indlist = [17 18 61 62 71 72];
indlist = [find([RhoStruct.motifpair.issame]==1 & strcmp({RhoStruct.motifpair.brainregion}, 'RA')) ...
    find([RhoStruct.motifpair.issame]==1 & strcmp({RhoStruct.motifpair.brainregion}, 'LMAN'))];
indlist = sort(indlist);

for j=1:length(indlist)
ind_motifpair = indlist(j);
bregionthis = RhoStruct.motifpair(ind_motifpair).brainregion;
rhoall = [];

% =============== 1) first do the first motif
motifthis = RhoStruct.motifpair(ind_motifpair).motif1;

% -- find all cases that have 1) this brain region, 2) one (and only one)
% of motifs is this single syl, 3) is not same type

inds = (strcmp({RhoStruct.motifpair.motif1}, motifthis) | strcmp({RhoStruct.motifpair.motif2}, motifthis)) ...
    & [RhoStruct.motifpair.issame]==0 & strcmp({RhoStruct.motifpair.brainregion}, bregionthis);
    
disp(['found N = ' num2str(sum(inds))]);

rhothis1 = [RhoStruct.motifpair(inds).rhovstime];

% =============== 2) then do second motif
motifthis = RhoStruct.motifpair(ind_motifpair).motif2;
% -- find all cases that have 1) this brain region, 2) one (and only one)
% of motifs is this single syl, 3) is not same type

inds = (strcmp({RhoStruct.motifpair.motif1}, motifthis) | strcmp({RhoStruct.motifpair.motif2}, motifthis)) ...
    & [RhoStruct.motifpair.issame]==0 & strcmp({RhoStruct.motifpair.brainregion}, bregionthis);
    
disp(['found N = ' num2str(sum(inds))]);

rhothis2 = [RhoStruct.motifpair(inds).rhovstime];

% ===== combine rho from 1st and 2nd motifs
rhoall = [rhothis1 rhothis2];

figure; hold on;
lt_subplot(1,2,1); hold on;
title(bregionthis);
plot(rhoall, '-r');
plot(mean(rhoall,2), '-r', 'LineWidth', 3);
plot(RhoStruct.motifpair(ind_motifpair).rhovstime, '-k', 'LineWidth', 3);
ylim([-1 1]);
lt_subplot(1,2,2); hold on;
% y = RhoStruct.motifpair(ind_motifpair).rhovstime-mean(rhoall,2);
y = RhoStruct.motifpair(ind_motifpair).rhovstime;
y = 2-(y+1);
plot(y, '-k', 'LineWidth', 3);
lt_plot_zeroline;
% ylim([-1 1]);
ylim([0 2]);
end

%% ========== PLOT ALL SAME AND DIFF PAIRS OVERLAYED
close all;

% ================ 1) separate plot for each brain region, overlay all rho
figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for bb =1 :length(BrainRegionList)
    bregionthis = BrainRegionList{bb};
    
    % --- all datapoints with this bregion
    indstmp = strcmp({RhoStruct.motifpair.brainregion}, bregionthis);
    
    RhoMat = [RhoStruct.motifpair(indstmp).rhovstime];
    IsSame = [RhoStruct.motifpair(indstmp).issame];
    
    
    % =============== 1) diff types (all trials)
    ymat = RhoMat(:, IsSame==0);
    plotcol = 'r';
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bregionthis ', diff']);
    ylabel('rho');
    xlabel('tbin');
    
    plot(1:size(ymat,1), ymat, '-', 'Color', plotcol);
    
    
    % =============== 2) same types (all trials)
    ymat = RhoMat(:, IsSame==1);
    plotcol = 'k';
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('same');
    
    plot(1:size(ymat,1), ymat, '-', 'Color', plotcol);
    
    
    % =========== 3) overlay mean and std of same and diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('same(k), diff (r)');
    
    % --- diff
    ymat = RhoMat(:, IsSame==0);
    plotcol = 'r';
    
    ymean = mean(ymat, 2);
    ystd = std(ymat, 0, 2);
    ysem = lt_sem(ymat');
        
%     shadedErrorBar(1:size(ymat,1), ymean, ystd, {'Color', plotcol}, 1);
    shadedErrorBar(1:size(ymat,1), ymean, ysem, {'Color', plotcol}, 1);
    
    % --- same
    ymat = RhoMat(:, IsSame==1);
    plotcol = 'k';
    
    
    ymean = mean(ymat, 2);
    ystd = std(ymat, 0, 2);
    ysem = lt_sem(ymat');
        
%     shadedErrorBar(1:size(ymat,1), ymean, ystd, {'Color', plotcol}, 1);
    shadedErrorBar(1:size(ymat,1), ymean, ysem, {'Color', plotcol}, 1);
    
    % -----
    lt_plot_zeroline;
end


%% =============== 2) SEPARATE PLOT FOR EACH PAIR --> COMPARE RA VS. LMAN
% ---- xcov parameters
maxlag = 0.06; % seconds
plotgaussfit = 0; % just for visualization, for all motifs
averageAcrossClasses = 1; % if 1, then averages over all motifpairs for a given single syl.

% ----
nummotifpairs = max([RhoStruct.motifpair.motifpairID]);
plotcollist = lt_make_plot_colors(length(BrainRegionList), 0,0);
figcount=1;
subplotrows=4;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
count = 1;
for j=1:nummotifpairs
    
    % --- check that have brain regions for this pair
    indtmp = find([RhoStruct.motifpair.motifpairID] == j);
    
    % --- sanity checks (make sure have all brain regions and no more than
    % that.
    assert(length(indtmp) == length(BrainRegionList),'asdfasd');
    assert(length(unique({RhoStruct.motifpair(indtmp).brainregion})) == length(indtmp), 'not all unique brain regions ..?');
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    
    issame = unique([RhoStruct.motifpair(indtmp).issame]); assert(length(issame) ==1, 'asd');
    if issame==1
        motif1 = RhoStruct.motifpair(indtmp(1)).motif1;
        motif2 = RhoStruct.motifpair(indtmp(1)).motif2;
       [~, singlesyl] =  lt_neural_QUICK_SylPairType(motif1, motif2);
    else
        singlesyl = '-';
    end
    
    if issame==1
        title(['[SAME]motifpair' num2str(j)]);
    else
        title(['[diff]motifpair' num2str(j)]);
    end
    
    for bb=1:length(BrainRegionList)
        bregionthis = BrainRegionList{bb};
        
        indtmp = strcmp({RhoStruct.motifpair.brainregion}, bregionthis) & ...
            [RhoStruct.motifpair.motifpairID] == j;
        assert(sum(indtmp) ==1,'asdf');
        
        y = RhoStruct.motifpair(indtmp).rhovstime;
        
        plot(y, '-', 'Color', plotcollist{bb});
        
    end
    ylim([-1 1]);
    lt_plot_zeroline;
    
    
    % ======== calculate cross correlation between traces for pairs of
    % brain regions
    for bb=1:length(BrainRegionList)
        bregion1 = BrainRegionList{bb};
        for bbb=bb+1:length(BrainRegionList)
            bregion2 = BrainRegionList{bbb};
            
            indtmp1 = strcmp({RhoStruct.motifpair.brainregion}, bregion1) & ...
                [RhoStruct.motifpair.motifpairID] == j;
            indtmp2 = strcmp({RhoStruct.motifpair.brainregion}, bregion2) & ...
                [RhoStruct.motifpair.motifpairID] == j;
            
            % ---- collect traces
            y1 = RhoStruct.motifpair(indtmp1).rhovstime;
            y2 = RhoStruct.motifpair(indtmp2).rhovstime;
            
            % --- get cross correlation, with brain regions in alphabetical
            % order
            [cc, lags] = xcov(y1, y2, floor(maxlag/binsize), 'coeff');
            [bregionout, inds] = sort({bregion1, bregion2});
            if inds(1)>inds(2)
                % -- then flip
                cc = flipud(cc);
            end
            
            % ====== OUTPUT
            % ----- 1) PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(bregionout);
            plot(lags, cc, '-k');
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            
            % -- overlay gaussian fit
            if plotgaussfit==1
            x = lags';
            y = cc;
            f = fit(x, y, 'gauss1');
            plot(f, x, y);

            end
            
            % ----- 2) SAVE
%             assert(length(RhoStruct.motifpairID) < j, 'overwriting data ...');
            RhoStruct.MotifpairBregionpair(count).bregion_inorder = bregionout;
            RhoStruct.MotifpairBregionpair(count).cc = cc;
            RhoStruct.MotifpairBregionpair(count).lags = lags;
            RhoStruct.MotifpairBregionpair(count).issame = issame;
            
            RhoStruct.MotifpairBregionpair(count).motifID = j;
            RhoStruct.MotifpairBregionpair(count).singlesyl = singlesyl;
            
            count = count+1;
            
        end
    end
end

%%  =================== PLOT ALL CC (separation in to same, diff)(also plot
% LMAN-LMAN, LMAN-RA, RA, RA)
BrainRegionPairs = {{'LMAN', 'LMAN'}, {'LMAN', 'RA'}, {'RA','RA'}};

bregionPair = {'LMAN', 'RA'};
AllIsSame = [];
AllCC = [];
AllLags = [];
AllSingleSyls = [];

for j=1:length(RhoStruct.MotifpairBregionpair)
    
    bregionthis = RhoStruct.MotifpairBregionpair(j).bregion_inorder;
    
    % -- is this bregion pair what I want?
    if ~all(strcmp(bregionPair, bregionthis))
        continue
    end
    
    % ---- collect this instance (bregion pair x motif pair)
    issame = RhoStruct.MotifpairBregionpair(j).issame;
    cc = RhoStruct.MotifpairBregionpair(j).cc;
    lags = RhoStruct.MotifpairBregionpair(j).lags;
    ssyl = RhoStruct.MotifpairBregionpair(j).singlesyl;
    
    AllIsSame = [AllIsSame issame];
    AllCC = [AllCC cc];
    AllSingleSyls = [AllSingleSyls ssyl];
        
%     AllLags = [AllLags lags'];   
end

% ============== COMBINE IF SAME SINGLE SYL
if averageAcrossClasses==1
   singlesyls = unique(AllSingleSyls);
   for syl = singlesyls
      
       if syl =='-'
           continue
       end
       
       numclasses = length(strfind(AllSingleSyls, syl));
       if numclasses ==1
           continue
       end
       
        disp(syl);
        
        indstoaverage = strfind(AllSingleSyls, syl);
        
        issame = unique(AllIsSame(indstoaverage))==1;
        assert(length(issame)==1, 'sdafsd');
        ytmp = AllCC(:, indstoaverage);
        
        ymean = mean(ytmp,2);
        
        % ----- remove these classes from dataset
    AllIsSame(indstoaverage) = [];
    AllCC(:, indstoaverage) = [];
    AllSingleSyls(indstoaverage) = [];
        
        % ----- add on the mean of these classes to dataset
    AllIsSame = [AllIsSame issame];
    AllCC = [AllCC ymean];
    AllSingleSyls = [AllSingleSyls syl];
   end
end

% ================ PLOT
lt_figure; hold on;
% ----- DIFF
lt_subplot(3,1,1); hold on;
plotcol = 'r';
plot(lags, AllCC(:, AllIsSame==0), '-', 'Color', plotcol);
plot(lags, mean(AllCC(:, AllIsSame==0), 2), 'LineWidth', 3, 'Color', plotcol);
lt_plot_zeroline;
lt_plot_zeroline_vert;
% ---- SAME
lt_subplot(3,1,2); hold on ;
plotcol = 'k';
plot(lags, AllCC(:, AllIsSame==1), '-', 'Color', plotcol);
plot(lags, mean(AllCC(:, AllIsSame==1), 2), 'LineWidth', 3, 'Color', plotcol);
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ---- COMBINED
lt_subplot(3,1,3); hold on;
% diff
plotcol = 'r';
plot(lags, mean(AllCC(:, AllIsSame==0), 2), 'LineWidth', 3, 'Color', plotcol);
% same
plotcol = 'k';
plot(lags, mean(AllCC(:, AllIsSame==1), 2), 'LineWidth', 3, 'Color', plotcol);

lt_plot_zeroline;
lt_plot_zeroline_vert;


% ==================== PLOT, guassian fits
lt_figure; hold on;
% % ----- DIFF
% lt_subplot(3,1,1); hold on;
% plotcol = 'r';
% plot(lags, AllCC(:, AllIsSame==0), '-', 'Color', plotcol);
% plot(lags, mean(AllCC(:, AllIsSame==0), 2), 'LineWidth', 3, 'Color', plotcol);
% lt_plot_zeroline;
% lt_plot_zeroline_vert;
% ---- SAME
lt_subplot(3,1,2); hold on ;
plotcol = 'k';

ccmat = AllCC(:, AllIsSame==1);
for j=1:size(ccmat,2)
   cc = ccmat(:,j);

    x = lags';
            y = cc;
            f = fit(x, y, 'gauss1');
            plot(f, x, y);
           legend(gca, 'off')
end
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ---- COMBINED
% lt_subplot(3,1,3); hold on;
% % diff
% plotcol = 'r';
% plot(lags, mean(AllCC(:, AllIsSame==0), 2), 'LineWidth', 3, 'Color', plotcol);
% % same
% plotcol = 'k';
% plot(lags, mean(AllCC(:, AllIsSame==1), 2), 'LineWidth', 3, 'Color', plotcol);
% 
% lt_plot_zeroline;
% lt_plot_zeroline_vert;
% 


        
%% ================ ASK HOW WELL A SINGLE DECODER CAN DECODE CONTEXT





    %% ===
    b = 3;
    indsneur = strcmp(NeurBrainRegions, 'LMAN');
    tmp = squeeze(ClassStruct.branch(b).NspkNeurTimeMotif(indsneur, 1, :));
        
        figure; hold on ; 
        plot(tmp(:,1), '-ok');
%         plot(frmat1(:,1), '-ob')
        plot(tmp(:,2), '-or');
        
        corr(tmp(:,1), tmp(:,2))

























