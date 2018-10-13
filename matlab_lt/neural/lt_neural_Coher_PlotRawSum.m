function lt_neural_Coher_PlotRawSum(COHSTRUCT, MOTIFSTATS_pop)


figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


numbirds = length(COHSTRUCT.bird);
for i=1:numbirds
    numexpt = length(COHSTRUCT.bird(i).experiment);
    birdname = SummaryStruct.birds(i).birdname;
    
    if ~any(strcmp(birdname, BirdsToPlot))
        continue
    end
    
    for ii=1:numexpt
        exptid = MOTIFSTATS_pop.birds(i).exptnum(ii).exptname;
        numsets = length(COHSTRUCT.bird(i).experiment(ii).setnum);
        
        for ss=1:numsets
            
            nummotifs = length(COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif);
            
            % ====== get list of neurons, bregions, and chans for this
            % dataset
            Neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{ss};
            Chanlist = [SummaryStruct.birds(i).neurons(Neurlist).channel];
            Bregionlist = {SummaryStruct.birds(i).neurons(Neurlist).NOTE_Location};
            
            for mm=1:nummotifs
                motifname = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss).motif(mm).regexpstr;
                
                CohCell = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Coh_ChpairByTrial;
                Chanpairs = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Chanpairs;
                tbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).t_relons;
                ffbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).ffbins;
                
                numchanpairs = size(Chanpairs, 1);
                
                for cc=1:numchanpairs
                    chansthis = Chanpairs(cc,:);
                    bregionsthis = {Bregionlist{Chanlist==chansthis(1)}, ...
                        Bregionlist{Chanlist==chansthis(2)}}; assert(length(bregionsthis)==2);
                    
                    % =================== COLLECT COHEROGRAMS ACROSS TRIALS
                    % --------- 1) COLLECT COHEROGRAMS FOR THIS PAIR
                    dim1 = size(CohCell{1},1);
                    dim2 = size(CohCell{1},2);
                    
                    % --------------- CONVERT CELLS TO A MATRIC
                    % use for loop because soemtime different size..
                    ntrials = size(CohCell,2);
                    cohmat = nan(length(tbins), length(ffbins), ntrials);
                    for tt=1:ntrials
                        if all(size(cohmat(:,:,tt))==size(CohCell{cc,tt}))
                            cohmat(:,:,tt) = CohCell{cc,tt};
                        else
                            disp('skipped! wrong size');
                        end
                    end
                    %                     cohmat = reshape(cell2mat(CohCell(cc, :)), dim1, dim2, []);
                    
                    % --------- 1) PLOT a single random trial
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title('single e.g.');
                    ylabel({[birdname '-' exptid], ...
                        ['set' num2str(ss) '-' motifname '-ch' num2str(chansthis) '[' bregionsthis{1} '-' bregionsthis{2} ']']});
                    indtmp = randi(size(cohmat,3), 1);
                    cohthis = cohmat(:,:,indtmp);
                    imagesc(tbins, ffbins, cohthis');
                    axis tight;
                    
                    % -------- 2) Plot trial average
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title('mean');
                    xlabel('time (rel alignment');
                    cohmean = nanmean(cohmat, 3);
                    imagesc(tbins, ffbins, cohmean');
                    axis tight;
                    
                    % -------- 3) line plots, separated by freq band
                    pcols = lt_make_plot_colors(length(ffbins), 1, [1 0 0]);
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title('each ff band');
                    
                    for k=1:length(ffbins)
                        
                        ffthis = ffbins(k);
                        cohthis = squeeze(cohmat(:,k,:));
                        
                        cmean = nanmean(cohthis,2);
                        csem = lt_sem(cohthis');
                        
                        plot(tbins, cmean, 'Color', pcols{k}, 'LineWidth', 1);
                        
                        % ---- plot text of ff
                        lt_plot_text(tbins(end), cmean(end), [num2str(ffthis)], ...
                            pcols{k}, 10);
                        
                        %                        fthis = ff(k);
                        %                        title(['ffbin:' num2str(fthis)]);
                        %                        hsplots = [hsplots hsplot];
                        %
                        %                        cohthis = squeeze(cohmat(:, i, :));
                        %                        ymean = mean(cohthis,2);
                        %                        ysem = lt_sem(cohthis');
                        %                        x = 1:length(ymean);
                        %                        lt_plot(x, ymean, {'Errors', ysem, 'Color', 'k'});
                    end
                    axis tight
                    ylim([0.2 1])
                    
                end
            end
            pause;
            close all;
        end
    end
end
