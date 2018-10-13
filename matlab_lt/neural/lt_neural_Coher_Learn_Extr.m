function OUTSTRUCT = lt_neural_Coher_Learn_Extr(SwitchStruct, SwitchCohStruct, plotON)
%% lt 10/9/18 - % THIS 1) plots and 2) collects to do stats


%%
% plotON=0; 

numbirds = length(SwitchCohStruct.bird);

%%
% =============== INITIATE

OUTSTRUCT.CohMat = {};
OUTSTRUCT.bnum = [];
OUTSTRUCT.enum = [];
OUTSTRUCT.switch = [];
OUTSTRUCT.motifnum = [];
OUTSTRUCT.issame = [];
OUTSTRUCT.istarg = [];
OUTSTRUCT.tvals = {};
OUTSTRUCT.indsbase = {};
OUTSTRUCT.indsWN = {};

% =============== RUN
for i=1:numbirds
    bthis = SwitchStruct.bird(i).birdname;
    numexpts = length(SwitchCohStruct.bird(i).exptnum);
    for ii=1:numexpts
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
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
            
            for mm=1:nummotifs
                motifthis = datthis.motifnum(mm).motifname;
                
                % ========= COLLECT MEAN COHERENCE PRE AND POST SWITCH
                cohmat = datthis.motifnum(mm).cohmat;
                tvals = datthis.motifnum(mm).tvals;
                indsbase = tvals>tstart & tvals<tswitch;
                indsWN = tvals>tswitch & tvals<tend;
                
                if isempty(cohmat)
                    continue
                end
                % --- what type of syl is this?
                istarg = 0;
                issame = nan;
                learndir = [];
                titlecol = 'k';
                if any(strcmp(targsyls, motifthis))
                    istarg =1;
                    learndir = targdirs(find(strcmp(targsyls, motifthis)));
                    titlecol = 'r';
                end
                
                if plotON==1
                    % ========= PLOT COHEROGRAMS
                    % -- baseline
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[BASE]' motifthis], 'Color', titlecol);
                    ylabel([bthis '-' ename '-sw' num2str(ss)]);
                    cohthis = cohmat(:,:,indsbase);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 1);
                    % -- WN
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[WN]'], 'Color', titlecol);
                    cohthis = cohmat(:,:,indsWN);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 1);
                    % =========== OVERLAY BASELINE AND WN FREQUENCY BANDS
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[base(dash), WN(solid)]'], 'Color',titlecol);
                    
                    cohthis = cohmat(:,:,indsbase);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 2, ':');
                    
                    cohthis = cohmat(:,:,indsWN);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 2, '-');
                    
                    % ========= NOTE DOWN SYLTYPE
                end
                
                % ######################## EXTRACT
                OUTSTRUCT.CohMat = [OUTSTRUCT.CohMat; cohmat];
                OUTSTRUCT.bnum = [OUTSTRUCT.bnum; i];
                OUTSTRUCT.enum = [OUTSTRUCT.enum; ii];
                OUTSTRUCT.switch = [OUTSTRUCT.switch; ss];
                OUTSTRUCT.motifnum = [OUTSTRUCT.motifnum; mm];
                OUTSTRUCT.issame = [OUTSTRUCT.issame; issame];
                OUTSTRUCT.istarg = [OUTSTRUCT.istarg; istarg];
                OUTSTRUCT.tvals = [OUTSTRUCT.tvals; tvals];
                OUTSTRUCT.indsbase = [OUTSTRUCT.indsbase; indsbase];
                OUTSTRUCT.indsWN = [OUTSTRUCT.indsWN; indsWN];
                
            end
        end
        if plotON==1
            pause; close all;
        end
    end
end

