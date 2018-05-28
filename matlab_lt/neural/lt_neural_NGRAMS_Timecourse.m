function OUTSTRUCT = lt_neural_NGRAMS_Timecourse(OUTSTRUCT, SummaryStruct, Params, ...
    PairTypesToCompare, nshufftmp, DoDecode, tcoursestyle)
%% lt 5/22/18 - Calculates FR diff, but in shifting window
% === to get timecourse of contextual modulation

onlyDoDesiredPairtypes =1; % if 0, does al paiors; if 1, then only desired
% list of pairtypes
%%
% nshufftmp = 10; % 50-100 is optimal, see DIAG below.
% NOTE: usually 100, but only matters for data in which
% want to zscore vs. self. for global z this just takes the
% first shuffle anyways.

% PairTypesToCompare = {'1  0  0', '1  1  1'};
% measure_to_recalc = 'absfrdiff'; % currently only 'absfrdiff' works

%%

% ------ convert window timing relative to data
motifpredur = Params.regexpr.motifpredur;
windowx = motifpredur+Params.window_prem;
strtype = Params.strtype;
% assert(strcmp(strtype, 'xaa')==1, 'have not coded for other string tyupes yet ...');

%% ===========
savedir = '/bluejay5/lucas/analyses/neural/NGRAMS';
savedir = [savedir '/' Params.dirname];

%% ========== go thru all birds and extract

% ============================= FIRST, move old version
% if strcmp(measure_to_recalc, 'absfrdiff')
%     OUTSTRUCT.All_AbsFRdiff_Zrelshuff_ORIG = OUTSTRUCT.All_AbsFRdiff_Zrelshuff;
%     OUTSTRUCT.All_AbsFRdiff_ORIG = OUTSTRUCT.All_AbsFRdiff;
%     OUTSTRUCT.All_AbsFRdiff_NEG_ORIG = OUTSTRUCT.All_AbsFRdiff_NEG;
%
%     OUTSTRUCT.All_AbsFRdiff_Zrelshuff = [];
%     OUTSTRUCT.All_AbsFRdiff = [];
%     OUTSTRUCT.All_AbsFRdiff_NEG = [];
%
%     OUTSTRUCT.All_N_ORIG = [];
%     OUTSTRUCT.All_N = [];
% end

OUTSTRUCT.AllTcourse_FRdiffDAT = cell(size(OUTSTRUCT.All_AbsFRdiff,1),1);
OUTSTRUCT.AllTcourse_FRdiffShuff = cell(size(OUTSTRUCT.All_AbsFRdiff,1),1);
OUTSTRUCT.AllTcourse_FRdiff_Z = cell(size(OUTSTRUCT.All_AbsFRdiff,1),1);
                

% ===================== SECOND, RUN
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    
    tmp = load([savedir '/bird' num2str(i) '.mat']);
    birdstruct = tmp.birdstruct;
    birdname = SummaryStruct.birds(i).birdname;
    
    % ============================== calculate things across all neurons
    numneurons = length(birdstruct.neuron);
    
    for nn=1:numneurons
        disp(' ############################################################### ');
        disp([birdname ', ' num2str(i) '-' num2str(nn)]);
        ngramlist = birdstruct.neuron(nn).ngramlist;
        
        
        % ========= go thru all pairs of ngrams.
        % ----------- all pairwise ngrams
        for j=1:length(ngramlist)
            if isempty(birdstruct.neuron(nn).ngramnum(j).DAT)
                continue
            end
            for jj=j+1:length(ngramlist)
                
                if isempty(birdstruct.neuron(nn).ngramnum(jj).DAT)
                    continue
                end
                
                % ===== skip this if has been thrown out because is b ad
                % syl
                indOutstruct = find(OUTSTRUCT.All_birdnum == i & OUTSTRUCT.All_neurnum==nn & ...
                    all(OUTSTRUCT.All_ngrampair_inorder == [j jj],2));
                if isempty(indOutstruct)
                    continue
                end
                assert(length(indOutstruct)<2, 'OUTSTRUCT and other dat arent matched');
                
                pairtype = OUTSTRUCT.All_diffsyl_PairType(indOutstruct);
                pairtype_logical = OUTSTRUCT.All_diffsyl_logical(indOutstruct, :);
                
                
                % ============ GET FRMAT across 2 motifs
                frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat;
                frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
                %                 N1 = size(frmat1,2);
                %                 N2 = size(frmat2, 2);
                %                 assert(N1>=Params.Nmin, 'asfsd');
                %                 assert(N2>=Params.Nmin, 'asfsd');
                
                % ==================== skip if is not correct pairtype
                if onlyDoDesiredPairtypes==1
                    ptypethis = OUTSTRUCT.PairTypesInOrder{pairtype};
                    if ~ismember(ptypethis, PairTypesToCompare)
                        continue
                    end
                end
                
              %% =========== Iterate over all time bins and get FR diff
                if strcmp(tcoursestyle, 'binned')
                % IN PROGRESS\
                
                % --- iterate over this:                
                [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth, FRmat, TrialInds, ...
                    DecodeDAT, DecodeShuff] = lt_neural_NGRAMS_QUICK_FRdiff(...
                    frmat1, frmat2, 1, windowx, Params.Nmin, nshufftmp, Params, ...
                    DoDecode);
                
                elseif strcmp(tcoursestyle, 'time')
                    
                 [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth, FRmat, TrialInds] ...
                     = lt_neural_NGRAMS_QUICK_FRdiffTcourse(frmat1, ...
                     frmat2, 1, nshufftmp, Params);
                    
                end
                
                OUTSTRUCT.AllTcourse_FRdiffDAT{indOutstruct} = FRdiffDAT;
                OUTSTRUCT.AllTcourse_FRdiffShuff{indOutstruct} = FRdiffShuff;
                OUTSTRUCT.AllTcourse_FRdiff_Z{indOutstruct} = FRdiff_Z;
                
                
                %% ======================= SANITY CHECK
                motif1 = birdstruct.neuron(nn).ngramlist{j};
                motif2 = birdstruct.neuron(nn).ngramlist{jj};
                % -- code: 111 means same at all 3 positions. 010 means
                % diff-same-diff, etc.
                diffsyl_logical = motif1~=motif2;
                assert(all(diffsyl_logical == pairtype_logical), 'not matched ...');
            end
            
            
            
        end
    end
end

