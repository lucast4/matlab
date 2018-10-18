function OUTSTRUCT = lt_neural_Coher_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls)
%% lt 10/9/18 - % THIS 1) plots and 2) collects to do stats

% removeBadSyls. leave at 1. this removes syls that follow target (i.e.
% acute effects. and also things that are "same-type" but actually a
% previous target..

%%
% plotON=0;
tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;
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

OUTSTRUCT.chanpair = [];
OUTSTRUCT.bregionpair = {};
OUTSTRUCT.motifname = {};
% =============== RUN
for i=1:numbirds
    bthis = SwitchStruct.bird(i).birdname;
    numexpts = length(SwitchCohStruct.bird(i).exptnum);
    for ii=1:numexpts
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
        foundatarg =0;
        for ss=1:numswitch
            if onlyfirstswitch==1
                if ss>1
                    continue
                end
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
            subplotcols=6;
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
                
                % ========= COLLECT MEAN COHERENCE PRE AND POST SWITCH
                cohmat = datthis.motifnum(mm).cohmat;
                % ------ take mean across all channels paiors
                if averagechanpairs ==1
                    cohmat = nanmean(cohmat, 4);
                end
                
                
                tvals = datthis.motifnum(mm).tvals;
                indsbase = tvals>tstart & tvals<tswitch;
                indsWN = tvals>tswitch & tvals<tend;
                
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
                if averagechanpairs==1
                    % then there is one series of cohmat, averaged over
                    % chanel p[airs
                    OUTSTRUCT.CohMat = [OUTSTRUCT.CohMat; cohmat];
                    OUTSTRUCT.bnum = [OUTSTRUCT.bnum; i];
                    OUTSTRUCT.enum = [OUTSTRUCT.enum; ii];
                    OUTSTRUCT.switch = [OUTSTRUCT.switch; ss];
                    OUTSTRUCT.motifnum = [OUTSTRUCT.motifnum; mm];
                    OUTSTRUCT.motifname = [OUTSTRUCT.motifname; motifthis]; 
                    OUTSTRUCT.issame = [OUTSTRUCT.issame; issame];
                    OUTSTRUCT.istarg = [OUTSTRUCT.istarg; istarg];
                    OUTSTRUCT.tvals = [OUTSTRUCT.tvals; tvals];
                    OUTSTRUCT.indsbase = [OUTSTRUCT.indsbase; indsbase];
                    OUTSTRUCT.indsWN = [OUTSTRUCT.indsWN; indsWN];
                else
                    % shoudl save one series for each channel pair
                    chanpairs = datthis.motifnum(mm).chanpair;
                    bregionpairs = datthis.motifnum(mm).bregionpair;
                    numpairs = size(chanpairs,1);
                    for k=1:numpairs
                        % -- save for this pair
                        OUTSTRUCT.CohMat = [OUTSTRUCT.CohMat; cohmat(:,:,:,k)];
                        OUTSTRUCT.bnum = [OUTSTRUCT.bnum; i];
                        OUTSTRUCT.enum = [OUTSTRUCT.enum; ii];
                        OUTSTRUCT.switch = [OUTSTRUCT.switch; ss];
                        OUTSTRUCT.motifnum = [OUTSTRUCT.motifnum; mm];
                        OUTSTRUCT.motifname = [OUTSTRUCT.motifname; motifthis]; 
                        OUTSTRUCT.issame = [OUTSTRUCT.issame; issame];
                        OUTSTRUCT.istarg = [OUTSTRUCT.istarg; istarg];
                        OUTSTRUCT.tvals = [OUTSTRUCT.tvals; tvals];
                        OUTSTRUCT.indsbase = [OUTSTRUCT.indsbase; indsbase];
                        OUTSTRUCT.indsWN = [OUTSTRUCT.indsWN; indsWN];
                        
                        OUTSTRUCT.chanpair = [OUTSTRUCT.chanpair; chanpairs(k,:)];
                        OUTSTRUCT.bregionpair = [OUTSTRUCT.bregionpair; bregionpairs{k}];
                    end
                    
                end
                
            end
        end
        if foundatarg==0 & numswitch>0
            % PROBLEM: this expt is missing motifs... see previous
            % extractio function, is likely eb cuase is throwing out
            % anything with nans.
            keyboard
        end
        if plotON==1
            pause; close all;
        end
    end
end


%% ======== 2) CALCULATE TRAIN MINUS BASE COHERENCE
cohdiff_all = {};
indsbase_epoch = {};
indsWN_epoch = {};
% --- takes 2nd half of training and base
for j=1:length(OUTSTRUCT.bnum)
    
    cohmat = OUTSTRUCT.CohMat{j};
    indsbase = find(OUTSTRUCT.indsbase{j});
    indsWN = find(OUTSTRUCT.indsWN{j});
    
    % -------------- take second half of WN inds
    indsbase = indsbase(round(length(indsbase)/2):end);
%     try
%     indsbase = indsbase(end-25:end);
%     catch err
%     end
    indsWN = indsWN(round(length(indsWN)/2):end);
    
    % -----
    cohmean_base = nanmean(cohmat(:,:, indsbase),3);
    cohmean_WN = nanmean(cohmat(:,:, indsWN),3);
    cohmean_diff = cohmean_WN - cohmean_base;
    
    if any(isnan(cohmean_diff(:)))
        keyboard
    end
    
    cohdiff_all = [cohdiff_all; cohmean_diff];
    indsbase_epoch = [indsbase_epoch; indsbase];
    indsWN_epoch = [indsWN_epoch; indsWN];
end
OUTSTRUCT.CohMean_WNminusBase = cohdiff_all;
OUTSTRUCT.indsbase_epoch = indsbase_epoch;
OUTSTRUCT.indsWN_epoch = indsWN_epoch;

