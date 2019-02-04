function [SwitchCohStruct, PARAMS] = lt_neural_LFP_PitchCorr(COHSTRUCT, SwitchCohStruct, SwitchStruct, ...
    PARAMS, twind, fwind, useWNtiming, WNprctile, prewind_relWN, interpol, ...
    onlyusegoodtargsyls)
%% lt 11/11/18 - scalar correaltion between coherence and pitch

if ~exist('useWNtiming', 'var')
    useWNtiming = 0;
end

if ~exist('onlyusegoodtargsyls', 'var')
    onlyusegoodtargsyls = 1;
end

%%

inds_f = PARAMS.ffbins>=fwind(1) & PARAMS.ffbins<=fwind(2);

%% EXTRACT SCALAR VALUE OF COHERENCE

numbirds = length(SwitchCohStruct.bird);
for i=1:numbirds
    bname = SwitchStruct.bird(i).birdname;
    numexpt = length(SwitchCohStruct.bird(i).exptnum);
    for ii=1:numexpt
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        for ss=1:numswitch
            nummotif = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum);
            
            if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum)
                continue
            end
            
            % === what does this translate into wrt to inds?
            if useWNtiming==0
                inds_t = PARAMS.tbins>=twind(1) & PARAMS.tbins<=twind(2);
            elseif useWNtiming==1
                % == get mean timing over all target syls
                targsyls = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningContingencies(1:2:end);
                prtilesall = [];
                for tt=1:length(targsyls)
                    sylthis = targsyls{tt};
                    
                    if onlyusegoodtargsyls==1
                        sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bname, ename, ss, sylthis);
                        if sylbad==1
                            continue
                        end
                    end
                    
                    indtmp = strcmp({SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum.motifname}, ...
                        sylthis);
                    if ~any(indtmp)
                        continue
                    end
                    assert(sum(indtmp)==1);
                    
                    dattmp = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(indtmp);
                    hittmp = dattmp.WNhittimes_min([dattmp.indsbase_epoch dattmp.indsWN_epoch]) - PARAMS.motif_predur;
                    hittmp = hittmp(~isnan(hittmp));
                    
                    prtilesall = [prtilesall prctile(hittmp, WNprctile)];
                end
                prtilesall = mean(prtilesall);
                twindtmp = prtilesall + prewind_relWN;
                inds_t = PARAMS.tbins>=twindtmp(1) & PARAMS.tbins<=twindtmp(2);
            end
            
            for mm=1:nummotif
                disp([num2str(i) '-' num2str(ii) '-' num2str(ss)]);
                datthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm);
                
                if isempty(datthis.lfpall)
                    continue
                end
                
                % ==================== extract coherene scalars for this dat
                fnamethis = [datthis.fileprefix '/Coh' datthis.filesuffix];
                pairstoget = datthis.chanpairstokeep;
                
                cohmat = lt_neural_LFP_QUICK_loadProcessDat(fnamethis, pairstoget);
                %                 cohmat =
                % ------------ COLLECT SCALAR
                if interpol==1
                    % NOT YET DONE - WILL HAVE TO INTERPOLATE SEPARATELY
                    % FOR EACH TRIAL, THAT IS A HASSLE...
                    %                     figure;
                    %                     surf(PARAMS.tbins, PARAMS.ffbins, cohmat)
                    
                else
                    cohscal = squeeze(nanmean(nanmean(cohmat(inds_t, inds_f, :,:),1),2));
                    
                end
                
                % ------------ SAVE OUTPUT
                SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).cohscalar = single(cohscal);
                
            end
        end
    end
end


%%  save params
PARAMS.cohscal_twind = twind;
PARAMS.cohscal_fwind = fwind;

