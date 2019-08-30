function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_daynight2(TrialStruct, ListOfExperiments, ...
    ParamsTrial, baseNormMethod, edgemethod, plotEachExptRaw)
%% day night, quick plotting, only experiments that were fully labeled.

%%

Nrends = 20; % at edges

%% figure out which experiments
if (0) % replaced by raw plotting below.
    % ===================== 1) PLOT EACH EXPERIMENT RAW
    numbirds = length(TrialStruct.birds);
    for i=1:numbirds
        numexpts = length(TrialStruct.birds(i).exptnum);
        bname = TrialStruct.birds(i).birdname;
        for ii=1:numexpts
            ename = TrialStruct.birds(i).exptnum(ii).exptname;
            isSDP = isfield(TrialStruct.birds(i).exptnum(ii).sylnum(1), 'INFO_SylDimensions');
            %     disp(isSDP)
            functmp = @(X)(strcmp(X{1}, bname) & strcmp(X{2}, ename));
            if ~any(cellfun(functmp, ListOfExperiments)) | isSDP==1
                % then skip
                continue
            end
            disp([bname '-' ename]);
            disp(isSDP)
            
            
            % ==================== PLOT
            %         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            %
            %         close all;
            ignoreDiffType=0;
            birdtoplot = bname;
            expttoplot = ename; % blank if don't care
            lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Raw(TrialStruct, ParamsTrial, ...
                ignoreDiffType, birdtoplot, expttoplot);
            %
            %
            %         % ================
            %         numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
            %     for ss =1:numsyls
            %
            %         sylthis = TrialStruct.birds(i).exptnum(ii).sylnum(ss)
            %
            %     end
            %
            
        end
    end
end
%% ################## COLLECT ACROSS EXPERIEMENTS
binedges = 6/24:(1/24):24/24; % in hours, start to end of day
numbirds = length(TrialStruct.birds);

All_FFvals = {};
All_Tvals = {};
All_bnum = [];
All_enum =[];
All_istarg = [];
All_issame = [];
All_sylnum = [];
All_dayedges = {};
All_daylearn_train = {};
All_nightlearn_train = {};
All_wnontime = [];
disp('NOTE: potentially throwing out syls that dont have baseline in same time bin');

for i=1:numbirds
    numexpts = length(TrialStruct.birds(i).exptnum);
    bname = TrialStruct.birds(i).birdname;
    for ii=1:numexpts
        ename = TrialStruct.birds(i).exptnum(ii).exptname;
        isSDP = isfield(TrialStruct.birds(i).exptnum(ii).sylnum(1), 'INFO_SylDimensions');
        %     disp(isSDP)
        functmp = @(X)(strcmp(X{1}, bname) & strcmp(X{2}, ename));
        if ~any(cellfun(functmp, ListOfExperiments)) | isSDP==1
            % then skip
            continue
        end
        disp([bname '-' ename]);
        disp(isSDP)
        
        % ================ STATS FOR THIS EXPT
        wnontime = TrialStruct.birds(i).exptnum(ii).WNontime;
        disp(['wnontime: ' num2str(wnontime)]);
        learndir = TrialStruct.birds(i).exptnum(ii).targlearndir;
        
        % ================
        numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        for ss =1:numsyls
            
            sylthis = TrialStruct.birds(i).exptnum(ii).sylnum(ss).syl;
            tvals = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
            ffvals = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
            istarg = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_istarget;
            issame = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_similar;
            
            % ============== SUBTRACT BASELINE CIRCADIAN
            % --- 1) get baseline circadian bins
            tbase = tvals(tvals<wnontime);
            ffbase = ffvals(tvals<wnontime);
            
            % ---
            if strcmp(baseNormMethod, 'circadian')
                tbase_withinday = tbase - floor(tbase);
                indsbins = discretize(tbase_withinday, binedges);
                ffbase_binned = grpstats(ffbase, indsbins, {'mean'});
                tbase_binned = unique(indsbins);
                
                % ---- 2) for all day, bin then subtract base
                tvals_withinday = tvals - floor(tvals);
                indsbins = discretize(tvals_withinday, binedges);
                assert(~any(isnan(indsbins)));
                [~, indstmp] = ismember(indsbins, tbase_binned);
                assert(sum(indstmp==0)/length(indstmp)<0.01, 'more than 1% files thrown out..');
                %             assert(all(tbase_binned(indstmp(1:end)) == indsbins));
                % OUTPUT:
                tvals = tvals(indstmp~=0); % thropwing out those without bins...
                ffvals = ffvals(indstmp~=0);
                indstmp = indstmp(indstmp~=0);
                % subtract base bins
                ffvals = ffvals - ffbase_binned(indstmp);
            elseif strcmp(baseNormMethod, 'onemean')
                % - one value for mean over all days
                ffvals = ffvals - mean(ffbase);
            end
            
            % ---- 3) for each day, flip if training is down
            if TrialStruct.FFalreadyFlippedLearnDir==0
            ffvals = ffvals*learndir;
            end
            
            % ---- 4) for each day, collect edge ff vals
            daylist = unique(floor(tvals));
            ffedges = nan(1, max(daylist)*2);
            for day = daylist'
                
                indstmp = floor(tvals)==day;
                ttmp = tvals(indstmp);
                ftmp = ffvals(indstmp);
                assert(all(diff(ttmp)>=0)); %  make sure is sorted
                
                if strcmp(edgemethod, 'edgemean')
                    ffedges(day*2-1) = mean(ftmp(1:Nrends));
                    ffedges(day*2) = mean(ftmp(end-Nrends+1:end));
                elseif strcmp(edgemethod, 'regression')
                    % -- fit regression
                    [b] =lt_regress(ftmp, ttmp, 0);
                    % -- collect FF at beginning and end of each day (from
                    % regression fit)
                    ff_fit = b(1) + b(2)*(ttmp);
                    ffedges(day*2-1) = ff_fit(1);
                    ffedges(day*2) = ff_fit(end);
                end
            end
            
            % --------- get overday and overnight, during train
            ffwn = ffedges(ceil(wnontime)*2-1:end); % only keep during WN on;
            ffday = ffwn(2:2:end) - ffwn(1:2:end);
            ffnight = ffwn(3:2:end) - ffwn(2:2:end-1);
            
            % ===== put back into struct;
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).ffedges =ffedges;
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).ffvals_minusbase =ffvals;
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).tvals_minusbase =tvals;
            
            
            
            % ============ OUTPUT
            All_FFvals = [All_FFvals; ffvals];
            All_Tvals = [All_Tvals; tvals];
            All_bnum = [All_bnum; i];
            All_enum =[All_enum; ii];
            All_istarg = [All_istarg; istarg];
            All_issame = [All_issame; issame];
            All_sylnum = [All_sylnum; ss];
            All_dayedges = [All_dayedges; ffedges];
            All_wnontime = [All_wnontime; wnontime];
            All_daylearn_train = [All_daylearn_train; ffday];
            All_nightlearn_train = [All_nightlearn_train; ffnight];
        end
        
    end
end


%% ================= PLOT EACH EXPERIEMNT, OPVERLAYING CALCULATED EDGE VALUES.
if plotEachExptRaw==1
    disp('ONLY WORKS WELL IF USE OVERALL BASE AS NORM...');
    maxbirds = max(All_bnum);
    maxexpts = max(All_enum);
    
    pcol = {};
    pcol{2,2} = [0.6 0.6 0.6]; %(targ, smae)
    pcol{1,2} = [0.3 0.3 0.8];
    pcol{1,1} = [0.8 0.3 0.3];
    
    for i=1:maxbirds
        for ii=1:maxexpts
            
            indsthis = All_bnum==i & All_enum==ii;
            if ~any(indsthis)
                continue
            end
            
            bname = TrialStruct.birds(i).birdname;
            ename = TrialStruct.birds(i).exptnum(ii).exptname;
            
            % ======= plot this experiment
            figcount=1;
            subplotrows=6;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            % ====== go thru all syls
            numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
            for ss=1:numsyls
                
                sylname = TrialStruct.birds(i).exptnum(ii).sylnum(ss).syl;
                t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).tvals_minusbase;
                ff_minusbase = TrialStruct.birds(i).exptnum(ii).sylnum(ss).ffvals_minusbase;
                ffedges = TrialStruct.birds(i).exptnum(ii).sylnum(ss).ffedges;
                ffedges = ffedges(~isnan(ffedges));
                tedges = reshape([unique(floor(t))'+0.2917; unique(floor(t))'+0.875], length(unique(floor(t)))*2, [])';
                assert(length(tedges)==length(ffedges));
                
                % is targ?
                istarg = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_istarget;
                issame = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_similar;
                pcthis = pcol{istarg+1, issame+1};
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title(sylname);
                if ss==1
                    ylabel([bname '-' ename]);
                end
                plot(t, ff_minusbase, 'x', 'Color', pcthis);
                lt_plot(tedges, ffedges);
                axis tight;
                lt_plot_zeroline;
            end
            
            % ====== plot summary for this expt
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            % hsplots = [hsplots hsplot];
            title('summary - actual analysis dat');
            % -- targ
            indtmp = All_bnum==i & All_enum==ii & All_istarg==1;
            ff = cell2mat(All_dayedges(indtmp));
            for j=1:length(ff)/2
                plot(j*2-1:j*2, ff(:, j*2-1:j*2)', '-xk');
            end
            
            % -- same
            indtmp = All_bnum==i & All_enum==ii & All_istarg==0 & ...
                All_issame==1;
            ff = cell2mat(All_dayedges(indtmp));
            for j=1:length(ff)/2
                plot(j*2-1:j*2, ff(:, j*2-1:j*2)', '-xb');
            end
            % -- diff
            indtmp = All_bnum==i & All_enum==ii & All_istarg==0 & ...
                All_issame==0;
            ff = cell2mat(All_dayedges(indtmp));
            for j=1:length(ff)/2
                plot(j*2-1:j*2, ff(:, j*2-1:j*2)', '-xr');
            end
            lt_plot_zeroline;
            
            % ----------------
            linkaxes(hsplots, 'xy');
        end
    end
end

%% ================ PLOT
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

indstoplot = lt_tools_grp2idx({All_bnum, All_enum});
for j = unique(indstoplot)'
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([]);
    
    % -- targ
    pcol = 'k';
    indsthis = indstoplot==j & All_istarg==1;
    y = cell2mat(All_dayedges(indsthis));
    wnontime = ceil(All_wnontime(indsthis));
    y = y(wnontime*2-1:end); % only keep during WN on;
    x = 1:length(y);
    plot(x, y, '-o', 'Color', pcol);
    plot(x(1:2:end), y(:, 1:2:end), 'ob');
    
    % --- same
    pcol = 'b';
    indsthis = indstoplot==j & All_istarg==0 & All_issame==1;
    y = cell2mat(All_dayedges(indsthis));
    wnontime = ceil(All_wnontime(indsthis));
    y = y(:, wnontime*2-1:end); % only keep during WN on;
    x = 1:size(y,2);
    plot(x, y', '-o', 'Color', pcol);
    
    % --- diff
    pcol = 'r';
    indsthis = indstoplot==j & All_istarg==0 & All_issame==0;
    y = cell2mat(All_dayedges(indsthis));
    wnontime = ceil(All_wnontime(indsthis));
    y = y(:, wnontime*2-1:end); % only keep during WN on;
    x = 1:size(y,2);
    plot(x, y', '-o', 'Color', pcol);
    
    
    % --
    lt_plot_zeroline;
    
end

lt_figure; hold on;
xlabel('targ, same, diff [day, night]');
ylabel('mean, each syl one val');

% ---- targ
x = [1 2];
pcol = 'k';
indsthis = All_istarg==1;
%
yday = cellfun(@nanmean, All_daylearn_train(indsthis));
ynight = cellfun(@nanmean, All_nightlearn_train(indsthis));
plot(x, [yday ynight], '-k');
lt_plot_bar(x, mean([yday ynight],1), {'Errors', lt_sem([yday ynight])})

% ---- same
x = [4 5];
pcol = 'b';
indsthis = All_istarg==0 & All_issame==1;
%
yday = cellfun(@nanmean, All_daylearn_train(indsthis));
ynight = cellfun(@nanmean, All_nightlearn_train(indsthis));
plot(x, [yday ynight], '-k');
lt_plot_bar(x, mean([yday ynight],1), {'Errors', lt_sem([yday ynight])})

% ---- diff
x = [7 8];
pcol = 'r';
indsthis = All_istarg==0 & All_issame==0;
%
yday = cellfun(@nanmean, All_daylearn_train(indsthis));
ynight = cellfun(@nanmean, All_nightlearn_train(indsthis));
plot(x, [yday ynight], '-k');
lt_plot_bar(x, mean([yday ynight],1), {'Errors', lt_sem([yday ynight])})









