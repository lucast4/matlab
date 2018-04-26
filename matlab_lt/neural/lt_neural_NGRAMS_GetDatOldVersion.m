 %% lt 4/24/18 - {OLD VERSION] to go from NGRAMSTRUCT to collected stats
 % NOW am not using that struct, since filesize tool arge to save. instead
 % am saving individual bird data then running lt_neural_NGRAMS_Compile
 % to extract stats.
 

% this code is stil usable.

 %%
    window_prem = [-0.03 0.03]; % relative to syl onset
Nshuffs = 1; % for negative control

% ------ convert window timing relative to data
motifpredur = NGRAMSTRUCT.params.Params.regexpr.motifpredur;
windowx = motifpredur+window_prem;
strtype = NGRAMSTRUCT.params.strtype;
assert(strcmp(strtype, 'xaa')==1, 'have not coded for other string tyupes yet ...');

% ==================
All_birdnum =[];
All_neurnum = [];
All_diffsyl_logical = [];
All_OneMinusRho = [];
All_OneMinusRho_NEG = [];

numbirds = length(NGRAMSTRUCT.bird);
for i=1:numbirds
    numneurons = length(NGRAMSTRUCT.bird(i).neuron);
    
    for nn=1:numneurons
        disp([num2str(i) '-' num2str(nn)])
        ngramlist = NGRAMSTRUCT.bird(i).neuron(nn).ngramlist;
        
        % ----------- all pairwise ngrams
        for j=1:length(ngramlist)
            if isempty(NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(j).DAT)
                continue
            end
            for jj=j+1:length(ngramlist)
                
                if isempty(NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(jj).DAT)
                    continue
                end
                
                % ##################################### RHO, FOR DATA
                % ================ COLLECT ALL TRIALS INTO ONE MATRIX
                segextract1 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(j).SegmentsExtract;
                segextract2 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(jj).SegmentsExtract;
                segextract1 = lt_neural_SmoothFR(segextract1);
                segextract2 = lt_neural_SmoothFR(segextract2);
                x = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(j).DAT.frx;
                
                frmat1 = [segextract1.FRsmooth_rate_CommonTrialDur];
                frmat2 = [segextract2.FRsmooth_rate_CommonTrialDur];
                if size(frmat1,1)~=size(frmat2,1)
                    indstmp = min([size(frmat1,1), size(frmat2,1)]);
                    frmat1 = frmat1(1:indstmp,:);
                    frmat2 = frmat2(1:indstmp,:);
                end
                FRmat = [frmat1 frmat2];
                TrialInds = {1:length(segextract1), ...
                    length(segextract1)+1:length(segextract1)+length(segextract2)};
                
                % --- pare down to just premotor windwo
                FRmat = FRmat(x>=windowx(1) & x<=windowx(2),:);                
                
                % ================= calculate correlation
                fr1 = mean(FRmat(:, TrialInds{1}),2);
                fr2 = mean(FRmat(:, TrialInds{2}),2);
                rhoDAT = corr(fr1, fr2);
                
                
                if (0) % OLD VERSION, cannot generalize to shuffle analyses.
                    % -------------- get correlation between premotor windows
                    fr1 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(j).DAT.frmean;
                    x1 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(j).DAT.frx;
                    fr1 = fr1(x1>=windowx(1) & x1<=windowx(2));
                    
                    fr2 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(jj).DAT.frmean;
                    x2 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(jj).DAT.frx;
                    fr2 = fr2(x2>=windowx(1) & x2<=windowx(2));
                    rho = corr(fr1, fr2);
                end
                
                % ################################## RHO, NEGATIVE CONTROL
                RhoShuff = [];
                nsamps = size(FRmat,2);
                for ss = 1:Nshuffs
                  indshuff = randperm(nsamps);
                  
                  FRmatSHUFF = FRmat(:,indshuff);
                  
                % ================= calculate correlation
                fr1 = mean(FRmatSHUFF(:, TrialInds{1}),2);
                fr2 = mean(FRmatSHUFF(:, TrialInds{2}),2);
                rho = corr(fr1, fr2);
                    
                RhoShuff = [RhoShuff rho];
                end
                rhoNEG = median(RhoShuff);
                
                % ------------- is this same pair or different pair?
                motif1 = NGRAMSTRUCT.bird(i).neuron(nn).ngramlist{j};
                motif2 = NGRAMSTRUCT.bird(i).neuron(nn).ngramlist{jj};
                % -- code: 111 means same at all 3 positions. 010 means
                % diff-same-diff, etc.
                diffsyl_logical = motif1~=motif2;
                
                
                
                % ################################## OUTPUT
                All_birdnum =[All_birdnum; i];
                All_neurnum = [All_neurnum; nn];
                All_diffsyl_logical = [All_diffsyl_logical; diffsyl_logical];
                All_OneMinusRho = [All_OneMinusRho; 1-rhoDAT];
                All_OneMinusRho_NEG = [All_OneMinusRho_NEG; 1-rhoNEG];
            end
        end
    end
end
