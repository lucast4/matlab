function lt_neural_v2_CTXT_Acoustic(analyfname, skipifdone)
%% lt 5/14/18 - extracts segmentsextract to pull out things that did not get before
% e.g. FF and other acoustic things (to add)

% skipifdone = 1 then skips indiviudal classes if already done.

%% lt 5/14/18 - confirm that decode is independent of FF correlation?
FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
% +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
collectWNhit = 0;

LearnKeepOnlyBase =1;


%% ========== load relevant dat
savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

load([savedir '/CLASSESv2_' analyfname '.mat']);
load([savedir '/SUMMARYv2_' analyfname '.mat']);
try
    load([savedir '/PARAMSv2_' analyfname '.mat']);
catch err
    
end

%% ========== get folder ready

savedirfinal = [savedir '/' analyfname '/FF'];
mkdir(savedirfinal);

%%
AllSizeMatch = []; % collect to determine whether size matches
AllBirdNeurBranchClass = [];
NumBirds = length(CLASSES.birds);
for i = 1:NumBirds
    Numneur = length(CLASSES.birds(i).neurons);
    
    for ii=1:Numneur
        
        
        % ------------------------- EXTRACT RAW DAT, THIS BIRD
        if isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
            % sam/mel data
            [SongDat, NeurDat, Params] = lt_neural_RASamMel_ExtractDat(SummaryStruct, i, ii);
            LearnKeepOnlyBase=0;
        else
            % my data
            [SongDat, NeurDat, Params] = lt_neural_ExtractDat2(SummaryStruct, i, ii);
        end
        
        % =============== GO THRU ALL BRANCHES
        Numbranch = length(CLASSES.birds(i).neurons(ii).branchnum);
        
        for bb=1:Numbranch
            
            
            % ====================== GO THRU ALL NEURONS
            
            %         DAT = ALLBRANCH.alignpos(apos).bird(i).branch(ii).neuron(nn);
            DAT = CLASSES.birds(i).neurons(ii).branchnum(bb);
            sylstr = DAT.regexprstr;
            
            Numclasses = length(DAT.SEGEXTRACT.classnum);
            disp(['bird' num2str(i) '_neur' num2str(ii) '_branch' num2str(bb)]);            
            for cc = 1:Numclasses
                
                % ================== if already done, then skip
                savefname = [savedirfinal '/bird' num2str(i) '_neur' num2str(ii) '_branch' num2str(bb) '_classnum' num2str(cc) '.mat'];
                if skipifdone==1
                if exist(savefname, 'file')~=0
                    % then already exists, skip analysis
                    disp(['[already done] skipping bird' num2str(i) '_neur' num2str(ii) '_branch' num2str(bb) '_classnum' num2str(cc)]);
                    continue
                end
                end
                   
                
                mclass = DAT.SEGEXTRACT.classnum(cc).regexpstr;
                
                % =========================== regexp to extract data
                % ---------------
                [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                    mclass, prms.motifpredur, prms.motifpostdur, prms.alignOnset, '', FFparams, ...
                    0, 1, collectWNhit, 0, LearnKeepOnlyBase, prms.preAndPostDurRelSameTimept);
                
                
                % --------------- test whether sample size is expected
                sizematch = length(SegmentsExtract) ==  length(DAT.SEGEXTRACT.classnum(cc).SegmentsExtract);
                AllSizeMatch = [AllSizeMatch; sizematch];
                AllBirdNeurBranchClass = [AllBirdNeurBranchClass; [i ii bb cc]];
                
                % ######################################## EXTRACT THINGS
                if isempty(SegmentsExtract)
                    continue
                end
                
                % ================== FF
                ff = [SegmentsExtract.FF_val];
                tvals = [SegmentsExtract.song_datenum];
                
                % ######################################## SAVE
%                 savefname = [savedirfinal '/bird' num2str(i) '_neur' num2str(ii) '_branch' num2str(bb) '_classnum' num2str(cc) '.mat'];
                
                % ------ save as matric [tvals, ff]; trials x 2
                t_ff = [tvals' ff'];
                save(savefname, 't_ff');
                
                
            end
        end
    end
end

% =================== SAVE VECTOR TELLING WHETHER DATA SIZE MATCHED
savefname = [savedirfinal '/AllSizeMatch.mat'];
save(savefname, 'AllSizeMatch');

savefname = [savedirfinal '/AllBirdNeurBranchClass.mat'];
save(savefname, 'AllBirdNeurBranchClass');




