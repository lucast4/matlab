function lt_neural_LFP_PlotEgRaw(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
    bnum, enum, swplot, motifnum, ntoplot, filt_low, filt_hi, SwitchCohStruct, ...
    SwitchStruct, savedir, saveON, plotCohScalExtremes, cohscaltmp, ...
    chanpairs, cohscal_allpairs)
%% lt 11/14/18 - plots output ofr lt_neural_LFP_PlotEgRaw_Extract
% NOTE: needs to run that function right before runnign this. run this as
% many times as desired for plotting examples. Can also save if desired.

if ~exist('plotCohScalExtremes', 'var')
    plotCohScalExtremes=0; 
end

% chanpairs, n x 2, n mnumber of pairs (array)
% cohscal_allpairs, cell nx1, each is 1xtrials.


i=bnum;
ii=enum;
mm = motifnum;
filt_fs = fs;

bname = SwitchStruct.bird(i).birdname;
ename = SwitchStruct.bird(i).exptnum(ii).exptname;

%%
datcohstruct = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm);

% -----
indsbase = datcohstruct.indsbase_epoch;
indsWN = datcohstruct.indsWN_epoch;
assert(size(datcohstruct.lfpall,2) == size(DatAll,1)) % n chans
assert(size(datcohstruct.lfpall,1) == size(DatAll{1},2)) % n trials

LFPdat = datcohstruct.lfpall;
t_lfp = datcohstruct.t_lfp;

if plotCohScalExtremes==1
    % 2 versions, either asks about all pairs and their coherence (top - acutlaly just 
    % takes the first channel pair) or older code where only care about one 
    % pair (doesn't work if multipel pairs exist for this experiment).
    % note, both cases uses bottom code.
    
    if exist('cohscal_allpairs', 'var')
        % NOTE: will look only at first pair. [i.e. to determine what
        % trials to take, will only take max and min for the first channel
        % pair]
        
        cohscaltmp = cohscal_allpairs{1}';

        xlabel_pairs = {};
        for j=1:size(chanpairs,1)
            xlabel_pairs = [xlabel_pairs [num2str(chanpairs(j,1)) '-' num2str(chanpairs(j,2))]];
        end
    end
    
    %  ========= PLOT LOW AND HIGH COHERENCE EXTREMES
    [~, indsort] = sort(cohscaltmp);
    
    indsort = indsort(ismember(indsort, [indsbase indsWN]));
    
    indslow = indsort(1:floor(ntoplot/2));
    indshi = indsort(end-ceil(ntoplot/2)+1:end);
    
    
%% ########################################### LOW
    figcount=1;
    subplotrows=length(chanlist_toget)+2;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    pcols_chans = lt_make_plot_colors(size(DatAll,1), 0, 0);
    for j=indslow'
        
        datfiltall = [];
        for cc = 1:length(chanlist_toget)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([bname '-' ename '-sw' num2str(swplot) '[LOW COH] trial' num2str(j) ',chan' num2str(chanlist_toget(cc)) ',' bregionlist{cc}]);
            lt_plot_annotation(1, ['coh=' num2str(cohscaltmp(j))], 'r');
            
            % =============== RAW DATA
            datmat = DatAll{cc};
            y = datmat(:, j);
            x = linspace(t_onoff{cc}(1), t_onoff{cc}(2), length(y));
            plot(x,y, 'Color', pcols_chans{cc});
            
            if ~exist('cohscal_allpairs', 'var')
            lt_plot_text(x(end), y(end), ['ch' num2str(chanlist_toget(cc))], pcols_chans{cc});
            end
            
            % ================= LFP
            datlfp = LFPdat{j,cc};
            %        plot(t_lfp, datlfp, 'Color', pcols_chans{cc}, 'LineWidth', 2);
            plot(t_lfp, datlfp, 'Color', 'k', 'LineWidth', 2);
            
            % =========== lfp (filtere)
            datfilt = lt_neural_filter(double(y), filt_fs, 0, filt_low, filt_hi);
            %        plot(x, datfilt, 'Color', pcols_chans{cc});
            plot(x, datfilt, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
            
            datfiltall = [datfiltall; datfilt'];
        end
        
        % ========= FILTER ALL AND OVERLAY
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        %     hsplots = [hsplots hsplot];
        for cc=1:size(DatAll,1)
            plot(x, datfiltall(cc,:), 'Color', pcols_chans{cc}, 'LineWidth', 2);
        end
        
        % ========= plot COHERENCE FOR ALL PAIRS FOR THIS TRIAL
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('coh scalar');
        xlabel('chan pair');
        cohthis = cellfun(@(x)x(j), cohscal_allpairs);
        lt_plot_bar(1:length(cohthis), cohthis);
        set(gca, 'XTick', 1:length(cohthis), 'XTickLabel', xlabel_pairs);
        ylim([0 1]);
        xlim([0 length(cohthis)+1]);
        
        
        axis tight;
        linkaxes(hsplots, 'xy');
    end
    
        
%% ########################################### LOW
    figcount=1;
    subplotrows=length(chanlist_toget)+2;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    pcols_chans = lt_make_plot_colors(size(DatAll,1), 0, 0);
    for j=indshi'
        
        datfiltall = [];
        for cc = 1:length(chanlist_toget)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([bname '-' ename '-sw' num2str(swplot) '[HI COH] trial' num2str(j) ',chan' num2str(chanlist_toget(cc)) ',' bregionlist{cc}]);
            lt_plot_annotation(1, ['coh=' num2str(cohscaltmp(j))], 'r');
            
            % =============== RAW DATA
            datmat = DatAll{cc};
            y = datmat(:, j);
            x = linspace(t_onoff{cc}(1), t_onoff{cc}(2), length(y));
            plot(x,y, 'Color', pcols_chans{cc});
            
            if ~exist('cohscal_allpairs', 'var')
            lt_plot_text(x(end), y(end), ['ch' num2str(chanlist_toget(cc))], pcols_chans{cc});
            end
            
            % ================= LFP
            datlfp = LFPdat{j,cc};
            %        plot(t_lfp, datlfp, 'Color', pcols_chans{cc}, 'LineWidth', 2);
            plot(t_lfp, datlfp, 'Color', 'k', 'LineWidth', 2);
            
            % =========== lfp (filtere)
            datfilt = lt_neural_filter(double(y), filt_fs, 0, filt_low, filt_hi);
            %        plot(x, datfilt, 'Color', pcols_chans{cc});
            plot(x, datfilt, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
            
            datfiltall = [datfiltall; datfilt'];
        end
        
        % ========= FILTER ALL AND OVERLAY
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        %     hsplots = [hsplots hsplot];
        for cc=1:size(DatAll,1)
            plot(x, datfiltall(cc,:), 'Color', pcols_chans{cc}, 'LineWidth', 2);
        end
        
        % ========= plot COHERENCE FOR ALL PAIRS FOR THIS TRIAL
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('coh scalar');
        xlabel('chan pair');
        cohthis = cellfun(@(x)x(j), cohscal_allpairs);
        lt_plot_bar(1:length(cohthis), cohthis);
        set(gca, 'XTick', 1:length(cohthis), 'XTickLabel', xlabel_pairs);
        ylim([0 1]);
        xlim([0 length(cohthis)+1]);
        
        
        axis tight;
        linkaxes(hsplots, 'xy');
    end
    
else

%% ########################################### BASE
figcount=1;
subplotrows=length(chanlist_toget)+1;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

indsbase_rand = indsbase(randperm(length(indsbase), min([ntoplot, length(indsbase)])));
pcols_chans = lt_make_plot_colors(size(DatAll,1), 0, 0);
for j=indsbase_rand
    
    datfiltall = [];
    for cc = 1:length(chanlist_toget)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title([bname '-' ename '-sw' num2str(swplot) '[BASE] trial' num2str(j) ',chan' num2str(chanlist_toget(cc)) ',' bregionlist{cc}]);
        
        % =============== RAW DATA
        datmat = DatAll{cc};
        y = datmat(:, j);
        x = linspace(t_onoff{cc}(1), t_onoff{cc}(2), length(y));
        plot(x,y, 'Color', pcols_chans{cc});
        lt_plot_text(x(end), y(end), ['ch' num2str(chanlist_toget(cc))], pcols_chans{cc});
        
        % ================= LFP
        datlfp = LFPdat{j,cc};
        %        plot(t_lfp, datlfp, 'Color', pcols_chans{cc}, 'LineWidth', 2);
        plot(t_lfp, datlfp, 'Color', 'k', 'LineWidth', 2);
        
        % =========== lfp (filtere)
        datfilt = lt_neural_filter(double(y), filt_fs, 0, filt_low, filt_hi);
        %        plot(x, datfilt, 'Color', pcols_chans{cc});
        plot(x, datfilt, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
        
        datfiltall = [datfiltall; datfilt'];
    end
    
    % ========= FILTER ALL AND OVERLAY
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     hsplots = [hsplots hsplot];
    for cc=1:size(DatAll,1)
        plot(x, datfiltall(cc,:), 'Color', pcols_chans{cc}, 'LineWidth', 2);
        
        
    end
    
    axis tight;
    linkaxes(hsplots, 'xy');
end



%% ########################################### WN
figcount=1;
subplotrows=length(chanlist_toget)+1;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

indsWN_rand = indsWN(randperm(length(indsWN), min([ntoplot, length(indsWN)])));
pcols_chans = lt_make_plot_colors(size(DatAll,1), 0, 0);
for j=indsWN_rand
    
    datfiltall = [];
    for cc = 1:length(chanlist_toget)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title([bname '-' ename '-sw' num2str(swplot) '[WN] trial' num2str(j) ',chan' num2str(chanlist_toget(cc)) ',' bregionlist{cc}]);
        
        % =============== RAW DATA
        datmat = DatAll{cc};
        y = datmat(:, j);
        x = linspace(t_onoff{cc}(1), t_onoff{cc}(2), length(y));
        plot(x,y, 'Color', pcols_chans{cc});
        lt_plot_text(x(end), y(end), ['ch' num2str(chanlist_toget(cc))], pcols_chans{cc});
        
        % ================= LFP
        datlfp = LFPdat{j,cc};
        %        plot(t_lfp, datlfp, 'Color', pcols_chans{cc}, 'LineWidth', 2);
        plot(t_lfp, datlfp, 'Color', 'k', 'LineWidth', 2);
        
        % =========== lfp (filtere)
        datfilt = lt_neural_filter(double(y), filt_fs, 0, filt_low, filt_hi);
        %        plot(x, datfilt, 'Color', pcols_chans{cc});
        plot(x, datfilt, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
        
        datfiltall = [datfiltall; datfilt'];
    end
    
    % ========= FILTER ALL AND OVERLAY
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     hsplots = [hsplots hsplot];
    for cc=1:size(DatAll,1)
        plot(x, datfiltall(cc,:), 'Color', pcols_chans{cc}, 'LineWidth', 2);
        
        
    end
    
    axis tight;
    linkaxes(hsplots, 'xy');
end

end
%% DISPLAY THE LIST OF CHANNEL PAIRS
for j=1:length(datcohstruct.bregionpair)
    disp([datcohstruct.bregionpair_originalorder{j} '  ===  ' num2str(datcohstruct.chanpair(j,:))]);
end

%% ======== save figures if needed
if saveON==1
    savedirthis = [savedir '/' bname '-' ename '-sw' num2str(swplot) '-motif' num2str(mm)];
    if ~exist(savedirthis)
        mkdir(savedirthis)
    else
        disp('OVERWRITING PREBVIOUS FIGURES!!');
    end
    cd(savedirthis);
    lt_save_all_figs;
end


