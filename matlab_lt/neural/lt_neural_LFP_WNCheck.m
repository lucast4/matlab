function lt_neural_LFP_WNCheck(SwitchCohStruct, MOTIFSTATS_pop, SwitchStruct, ...
    PARAMS)
%% lt 12/10/18 - across all switches and all syls, checks 1) fraction and 2)timing of WN hits

%% ######################### FRACTIONS WN HITS

figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:length(SwitchCohStruct.bird)
    bname = SwitchStruct.bird(i).birdname;
    for ii=1:length(SwitchCohStruct.bird(i).exptnum)
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        for ss=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
            motifnames ={};
            istarg = [];
            frachit = [];
            
            for mm=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum)
                
                if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).motifname)
                    continue
                end
                motifthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).motifname;
                neurset = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).neursetused;
                segextract = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(neurset).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                t = any(strcmp(motifthis, SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs(1:2:end)));
                
                trialsthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsbase_epoch(1):...
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsWN_epoch(end);
                
                fh = sum([segextract(trialsthis).hit_WN])./length(segextract(trialsthis));
                
                % === COLLECT
                motifnames = [motifnames; motifthis];
                frachit = [frachit; fh];
                istarg = [istarg; t];
            end
            
            if isempty(frachit)
                continue
            end
            
            % ======= plot
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([bname '-' ename '-s' num2str(ss)]);
            x = 1:length(frachit);
            lt_plot_bar(x, frachit);
            lt_plot(x(istarg==1), 0.9, {'Color', 'r'});
            lt_plot_bar(x(istarg==1), frachit(istarg==1), {'Color', 'r'});
            set(gca, 'XTick', x);
            set(gca, 'XTickLabel', motifnames);
            rotateXLabels(gca, 90);
            ylim([0 1]);
        end
    end
end


%% ######################### TIMES OF WN HITS

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:length(SwitchCohStruct.bird)
    bname = SwitchStruct.bird(i).birdname;
    for ii=1:length(SwitchCohStruct.bird(i).exptnum)
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        for ss=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
            motifnames ={};
            istarg = [];
            frachit = [];
            
            for mm=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum)
                
               if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).motifname)
                    continue
                end

                motifthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).motifname;
                neurset = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).neursetused;
                segextract = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(neurset).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                t = any(strcmp(motifthis, SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs(1:2:end)));
                                
                trialsthis = min(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsbase_epoch):...
                    max(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsWN_epoch);
                
                fh = sum([segextract(trialsthis).hit_WN])./length(segextract(trialsthis));

                if fh<0.05 & t==0
                    continue
                end
               
                wnonsets = [segextract(trialsthis).WNonset_sec];
                % -- make relative to syl onset
                wnonsets = wnonsets - PARAMS.motif_predur;
                
                if t==0
                    % then is not tareget, ignore
                    continue
                end
                
                % ===== PLOT WN ONSET DISTYRIBUTISON
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                if t==1
                title({[bname '-' ename '-s' num2str(ss)], [motifthis '[TARG]']});
                else
                    title({[bname '-' ename '-s' num2str(ss)], [motifthis '[NONTARG]']});
                end
                %             xcenters = -PARAMS.motif_predur:0.01:
                xcenters = -0.02:0.005:0.05;
                lt_plot_histogram(wnonsets, xcenters, 1, 0, '', 0, 'r');
                axis tight
                XLIM = xlim;
                %             xlim([-PARAMS.motif_predur XLIM(2)]);
                lt_plot_zeroline_vert;
                xlabel('onset time (rel suyl onset)');
            end
        end
    end
end