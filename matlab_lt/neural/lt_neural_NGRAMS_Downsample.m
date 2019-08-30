function [OUTSTRUCT_subsamp, downfactorlist] = lt_neural_NGRAMS_Downsample(OUTSTRUCT, SummaryStruct, Params, ...
    measure_to_recalc, downfactorlist)

nshufftmp = 2; % for getting negative control.
% downfactorlist = 0.1:0.1:1; % fractin of samples between N and Nmin.

assert(strcmp(measure_to_recalc, 'absfrdiff'), 'have not coded other measures...');

%% lt 5/4/18 - multiple downsamples to try to match RA data to LMAN
% since lMAN more noisy, so need to have lower RA sample to match?


%% ===========

savedir = '/bluejay5/lucas/analyses/neural/NGRAMS';
savedir = [savedir '/' Params.dirname];

% ---------- premtoor iwnodw tmiing
motifpredur = Params.regexpr.motifpredur;
windowx = motifpredur+Params.window_prem;


%% ================= go thru all pairtypes and downsample

% prepare out, will have 1st d indices matched to OUTSTRUCT.
OUTSTRUCT_subsamp.FRdiffDAT = nan(length(OUTSTRUCT.All_AbsFRdiff), length(downfactorlist));
OUTSTRUCT_subsamp.FRdiffShuff = nan(length(OUTSTRUCT.All_AbsFRdiff), length(downfactorlist));
OUTSTRUCT_subsamp.FRdiff_Z = nan(length(OUTSTRUCT.All_AbsFRdiff), length(downfactorlist));
OUTSTRUCT_subsamp.Nboth = nan(length(OUTSTRUCT.All_AbsFRdiff), 2, length(downfactorlist));


numbirds = length(SummaryStruct.birds);
counter = 1;
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
        
        
        % ----------- all pairwise ngrams
        for j=1:length(ngramlist)
            if isempty(birdstruct.neuron(nn).ngramnum(j).DAT)
                continue
            end
            for jj=j+1:length(ngramlist)
                
                if isempty(birdstruct.neuron(nn).ngramnum(jj).DAT)
                    continue
                end
                
                % ================ skip this if have already removed
                % because is a bad syl
                indtmp = find(OUTSTRUCT.All_birdnum == i & OUTSTRUCT.All_neurnum==nn & ...
                    all(OUTSTRUCT.All_ngrampair_inorder == [j jj],2));
                if isempty(indtmp)
                    % then this has been removed
                    continue
                end
                
                assert(indtmp == counter, 'saved dat and outstruct not matched ..')
                assert(length(indtmp)<2, 'OUTSTRUCT and other dat arent matched');
                
                
                
                % =================== CALCULATE STATS AS FUNCTION OF
                % DOWNSAMPLED SIZE
                % ------------ 1) GET FRMAT across 2 motifs
                frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat;
                frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
                
                % -------------- 2) get FR differences, for different
                % downsample factors
                FRdiffstruct = struct;
                for k=1:length(downfactorlist)
                    downfactor = downfactorlist(k); % to downsample, while still at least equal to Nmin
                    [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth] = lt_neural_NGRAMS_QUICK_FRdiff(...
                        frmat1, frmat2, downfactor, windowx, Params.Nmin, nshufftmp, ...
                        Params, Params.dodecode);
                    
                    % ========= save
                    FRdiffstruct(k).FRdiffDAT = FRdiffDAT;
                    FRdiffstruct(k).FRdiffShuff = FRdiffShuff;
                    FRdiffstruct(k).FRdiff_Z = FRdiff_Z;
                    FRdiffstruct(k).Nboth= Nboth;
                end
                
                
                % ============ plot function of fr diff vs. sample size
                if (0)
                    lt_figure; hold on;
                    title('data(k), neg (r)');
                    % --- dat
                    x = mean([FRdiffstruct.Nboth],1);
                    y = [FRdiffstruct.FRdiffDAT];
                    plot(x, y, '-ok');
                    
                    % --- neg
                    y = [FRdiffstruct.FRdiffShuff];
                    plot(x, y', '-r');
                    
                end
                
                
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 OUTPUT, SAVE ALL
                % DEOWNSAMPEL ITERATIONS
                OUTSTRUCT_subsamp.FRdiffDAT(counter, :) = [FRdiffstruct.FRdiffDAT];
                tmp = [FRdiffstruct.FRdiffShuff]; % just take one shuffle iteration
                OUTSTRUCT_subsamp.FRdiffShuff(counter, :) = tmp(1,:);
                OUTSTRUCT_subsamp.FRdiff_Z(counter, :) = [FRdiffstruct.FRdiff_Z];
                OUTSTRUCT_subsamp.Nboth(counter, :, :) = [FRdiffstruct.Nboth];
                
                % ===================
                counter = counter +1; % this should correspond to current position in OUTSTRUCT
                
            end
        end
    end
end

assert(counter == length(OUTSTRUCT.All_AbsFRdiff)+1, 'indices dont match...');








