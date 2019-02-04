function lt_neural_Coher_PlotChanMotif(SwitchStruct, SwitchCohStruct, OUTSTRUCT, ...
    bird, expt, sw, averagechanpairs, PARAMS, removeBadSyls, ffbinsedges)
%% lt 10/9/18 - % THIS 1) plots

%%
% plotON=0;
try
    tbins = PARAMS.tbins_BeforeAlignWN;
catch err
tbins = PARAMS.tbins;
end
ffbins = PARAMS.ffbins;
numbirds = length(SwitchCohStruct.bird);
%%
% =============== RUN
for i=1:numbirds
    bthis = SwitchStruct.bird(i).birdname;
    numexpts = length(SwitchCohStruct.bird(i).exptnum);
    
    for ii=1:numexpts
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
        foundatarg =0;
        for ss=1:numswitch
            
            
            if ~(strcmp(bird, bthis) & strcmp(expt, ename) & ss==sw)
                continue
            end
            
            nummotifs = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum);
            datthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss);
            if isempty(datthis.motifnum)
                continue
            end
            % =========== for this switch, get edge times and switch time
            tstart = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_previous;
            tswitch = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
            tend = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_next;
            targsyls = {SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{1:2:end}};
            targdirs = [SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{2:2:end}];
            motiflist = {datthis.motifnum.motifname};
            motiflist = motiflist(~cellfun(@isempty, motiflist));
            sylssame = lt_neural_QUICK_ExtractSameType(motiflist, targsyls);
            
            % =========== one figure for each switch
            figcount=1;
            subplotrows=3;
            subplotcols=8;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            if isempty(targsyls)
                keyboard
            end
            
            %             disp(['TARG SYLS: ' num2str(targsyls)]);
            for mm=1:nummotifs
                motifthis = datthis.motifnum(mm).motifname;
                
                
                if removeBadSyls==1
                    sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bthis, ename, ss, motifthis);
                    if sylbad==1
                        continue
                    end
                end
                
                
                % ============= FROM OUTSTRUCT
                indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss ...
                    & OUTSTRUCT.motifnum==mm);
                if isempty(indsthis)
                    continue
                end
                indsbase = OUTSTRUCT.indsbase_epoch{indsthis(1)};
                indsWN = OUTSTRUCT.indsWN_epoch{indsthis(1)};
                
                % ========= COLLECT MEAN COHERENCE PRE AND POST SWITCH
                if isfield(datthis.motifnum(mm), 'cohmat')
                    cohmat = datthis.motifnum(mm).cohmat;
                else
                    cohmat = lt_neural_LFP_QUICK_loadProcessDat([datthis.motifnum(mm).fileprefix '/Coh' datthis.motifnum(mm).filesuffix], ...
                        datthis.motifnum(mm).chanpairstokeep, 0);
                end
                chanpairs = datthis.motifnum(mm).chanpair;
                
                % ------ take mean across all channels paiors
                if averagechanpairs ==1
                    cohmat = nanmean(cohmat, 4);
                end
                
                
                tvals = datthis.motifnum(mm).tvals;
                
%                 indsbase = tvals>tstart & tvals<tswitch;
%                 indsWN = tvals>tswitch & tvals<tend;
%                 
                if isempty(cohmat)
                    continue
                end
                
                % --- TARG? what type of syl is this?
                if any(strcmp(targsyls, motifthis))
                    istarg =1;
                    learndir = targdirs(find(strcmp(targsyls, motifthis)));
                    titlecol = 'r';
                else
                    istarg = 0;
                    learndir = [];
                    titlecol = 'k';
                end
                
                % --- SAME?
                issame = any(strcmp(sylssame, motifthis));
                
                if istarg==1
                    foundatarg=1;
                end
                
                % ========= PLOT COHEROGRAMS
                % ---------------- 1 plot for each chan
                numchans = size(cohmat,4);
                for cc = 1:numchans
                    if averagechanpairs==0
                    chanpair = chanpairs(cc,:);
                    else
                        chanpair = 'average';
                    end
                    % -- baseline
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[BASE]' motifthis], 'Color', titlecol);
                    ylabel({[bthis '-' ename '-sw' num2str(ss)], ['chpair:' num2str(chanpair)]});
                    cohthis = cohmat(:,:,indsbase, cc);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 1, '', [], 0, 0, ffbinsedges);
%                     colorbar('East')
lt_plot_colormap('pval');
                    % -- WN
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[WN]'], 'Color', titlecol);
                    cohthis = cohmat(:,:,indsWN, cc);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 1, '', [], 0, 0, ffbinsedges);
lt_plot_colormap('pval');
%                     colorbar('East')
                    % == WN minus base
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[WN-base]'], 'Color', titlecol);
                    cohthis = nanmean(cohmat(:,:,indsWN, cc),3) - nanmean(cohmat(:,:,indsbase, cc),3);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 1, '', [-0.15 0.15], 0, 0, ffbinsedges);
%                     colorbar('East')
                    % =========== OVERLAY BASELINE AND WN FREQUENCY BANDS
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[base(dash), WN(solid)]'], 'Color',titlecol);
                    
                    cohthis = cohmat(:,:,indsbase, cc);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 2, ':', [], 0, 0, ffbinsedges);
                    
                    cohthis = cohmat(:,:,indsWN, cc);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 2, '-', [], 0, 0, ffbinsedges);
%                     lt_neural_Coher_Plot(cohthis, tbins, ffbins, 2, '-');
                    
                    % ========= NOTE DOWN SYLTYPE
                end
            end
        end
    end
end
