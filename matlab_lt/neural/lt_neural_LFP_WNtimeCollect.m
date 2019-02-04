function lt_neural_LFP_WNtimeCollect(SwitchCohStruct, MOTIFSTATS_pop, SwitchStruct, ...
    PARAMS)
%% LT 1/4/19 - extract time of WN hits for each moptif, save in outstruct.
% NOTE: ABANDONED !!! SINCE I REALIZED THAT THIS IS ALREADY EXTRACTED IN
% SWITCHCOHSTRUCT

%% ######################### TIMES OF WN HITS

for i=1:length(SwitchCohStruct.bird)
    %     bname = SwitchStruct.bird(i).birdname;
    for ii=1:length(SwitchCohStruct.bird(i).exptnum)
        %         ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        for ss=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
            %             motifnames ={};
            %             istarg = [];
            %             frachit = [];
            
            for mm=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum)
                
                if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).motifname)
                    continue
                end
                
                motifthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).motifname;
                neurset = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).neursetused;
                
                segextract = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(neurset).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
%                 t = any(strcmp(motifthis, SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs(1:2:end)));
                
                segextract.WNonset_sec

                % STOPPED HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!
                trialsthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsbase_epoch(1):...
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsWN_epoch(end);
                
                fh = sum([segextract(trialsthis).hit_WN])./length(segextract(trialsthis));
                
                wnonsets = [segextract(trialsthis).WNonset_sec];
                % -- make relative to syl onset
                wnonsets = wnonsets - PARAMS.motif_predur;
                
            end
        end
    end
end