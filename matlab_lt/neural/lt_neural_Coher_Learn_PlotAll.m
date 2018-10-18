for i=2:numbirds
    bthis = SwitchStruct.bird(i).birdname;
    numexpts = length(SwitchCohStruct.bird(i).exptnum);
    for ii=1:numexpts
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
        foundatarg =0;
        for ss=1:numswitch
            nummotifs = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum);
            datthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss);
            
            % =========== for this switch, get edge times and switch time
            tstart = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_previous;
            tswitch = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
            tend = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_next;
            targsyls = {SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{1:2:end}};
            targdirs = [SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{2:2:end}];
            
            % =========== one figure for each switch
            figcount=1;
            subplotrows=3;
            subplotcols=6;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            if isempty(targsyls)
                keyboard
            end
            
            for mm=1:nummotifs
                motifthis = datthis.motifnum(mm).motifname;
                
                % ========= COLLECT MEAN COHERENCE PRE AND POST SWITCH
                cohmat = datthis.motifnum(mm).cohmat;
                if isempty(cohmat)
                    continue
                end
                
                % --------- extract epoch bins.
                indsbase = datthis.motifnum(mm).indsbase_epoch;
                indsWN = datthis.motifnum(mm).indsWN_epoch;
%                 tvals = datthis.motifnum(mm).tvals;
%                 indsbase = find(tvals>tstart & tvals<tswitch);
%                 indsWN = tvals>tswitch & tvals<tend;
                
                % --- what type of syl is this?
                if any(strcmp(targsyls, motifthis))
                    istarg =1;
                    learndir = targdirs(find(strcmp(targsyls, motifthis)));
                    titlecol = 'r';
                else
                    istarg = 0;
                    issame = nan;
                    learndir = [];
                    titlecol = 'k';
                end
                
                numchanpairs = size(cohmat,4);
                
                % ==================== PLOT EACH CHANNEL PAIR
                for cc = 1:numchanpairs
                    
                    pairthis = datthis.motifnum(mm).chanpair(cc,:);
                                        
                    % ========= PLOT COHEROGRAMS
                    
                    % -- baseline
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[BASE]' motifthis], 'Color', titlecol);
                    ylabel({[bthis '-' ename '-sw' num2str(ss)], ['pair:' num2str(pairthis)]});
                    cohthis = cohmat(:,:,indsbase, cc);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 1);
                    
                    % -- WN
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[WN]'], 'Color', titlecol);
                    cohthis = cohmat(:,:,indsWN, cc);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 1);
                    
                    % =========== OVERLAY BASELINE AND WN FREQUENCY BANDS
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[base(dash), WN(solid)]'], 'Color',titlecol);
                    
                    cohthis = cohmat(:,:,indsbase, cc);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 2, ':');
                    
                    cohthis = cohmat(:,:,indsWN, cc);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 2, '-');
                    
                end
            end
            
            pause; close all;
            
        end
    end
end
