function DatStructByBranch = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH, ...
    analyfname, tbin, useDprime, timewindow, prms)
% timewindow = [-0.03 0.03]; % relative to syllable onset (for getting
% other metric of fr difference, one scalar for premotor window)

if ~exist('timewindow', 'var')
    timewindow = [];
end

%%  lt 11/14/17 - if also want to extract premotor window decoding (shuffle), then need:



%% lt 11/29/17 - organizes data by Unique branches (defined by regexpr str)
% OUTPUTS A NEW STRUCTURE


%%
apos = 1; % assume is new analysis version, this always 1.

if ~exist('analyfname', 'var')
    analyfname = '';
end

if ~exist('tbin', 'var')
    tbin = [];
end

if ~exist('useDprime', 'var')
    useDprime=0;
end

%% extract all branches first

ALLBRANCH = lt_neural_v2_CTXT_GetListOfBranches(ALLBRANCH);


%% PRODUCE NEW STRUCTURE ORGANIZED BY UNIQUE REGEXP STR

numbirds = length(ALLBRANCH.alignpos(apos).bird);
DatStructByBranch = struct;

for i=1:numbirds
    
    ListOfBranches = ALLBRANCH.alignpos(apos).bird(i).ListOfBranches;
    
    
    % ################################### initiate structure
    for kk=1:length(ListOfBranches)
        DatStructByBranch.bird(i).branchID(kk).regexpstr = '';
        DatStructByBranch.bird(i).branchID(kk).DAT.xdecode = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.ydecode = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.ydecode_neg = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.ydecode_pos = [];
        
        DatStructByBranch.bird(i).branchID(kk).DAT.y_decode_z = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.y_decode_negshuff_mean = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.y_decode_negshuff_std = [];
        
        %         DatStructByBranch.bird(i).branchID(kk).DAT.dprime = [];
        %         DatStructByBranch.bird(i).branchID(kk).DAT.dprime_neg = [];
        %
        DatStructByBranch.bird(i).branchID(kk).DAT.sylcontour_mean = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.sylcontour_x = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.brainregion = {};
        DatStructByBranch.bird(i).branchID(kk).DAT.PREMOTORDECODE_pval = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.PREMOTORDECODE_struct = [];
        
        DatStructByBranch.bird(i).branchID(kk).DAT.IndsOriginal_apos = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.IndsOriginal_bird = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.IndsOriginal_branch = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.IndsOriginal_neuron = [];
        
        DatStructByBranch.bird(i).branchID(kk).DAT.FRcorr_dat = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.FRcorr_neg = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.FRcorr_pos = [];
        
        DatStructByBranch.bird(i).branchID(kk).DAT.FFstruct = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.frdat = [];
    end
    
    
    % ################################### COLLECT DATA ND PUT INTO STRUCTURE
    numbranches = length(ALLBRANCH.alignpos(apos).bird(i).branch);
    for bb =1:numbranches
        
        numneurons = length(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron);
        for nn=1:numneurons
            
            % ================ collect data
            x_decode = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).xtimes;
            y_decode = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals;
            thisbranch = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).prms_regexpstr;
            
            y_decode_neg = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals_neg;
            y_decode_pos = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals_pos;
            
            sylcontour_mean = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).sylcontours_mean;
            sylcontour_x = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).sylcontours_x;
            
            % -------------------- keep data?
            if isempty(x_decode)
                continue
            end
            
            
            % ============ decode, z-scored to shuffle distribution
            if isfield(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn), 'yvals_z')
                y_decode_z = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals_z;
                y_decode_negshuff_mean = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals_neg_shuffmean;
                y_decode_negshuff_std = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals_neg_shuffstd;
            else
                y_decode_z = nan(size(y_decode));
                y_decode_negshuff_mean =  nan(size(y_decode));
                y_decode_negshuff_std =  nan(size(y_decode));
            end
            
            
            % ============= USE DPRIME? --- THEN REPLACE DATA
            if useDprime==1
                
                disp(size(y_decode));
                
                dprime = mean(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).DprimeAllPairwise,2); % mean dprime across branches
                dprime_neg = mean(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).DprimeAllPairwise_Neg, 2);
                dprime_pos = mean(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).DprimeAllPairwise_Pos, 2);
                
                disp(size(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).DprimeAllPairwise));
                
                x_decode = sylcontour_x;
                y_decode = dprime';
                y_decode_neg = dprime_neg';
                y_decode_pos = dprime_pos';
                
            end
            
            % ------------------- premotor window shuffles
            if ~isempty(analyfname)
                decodestruct = lt_neural_v2_CTXT_BRANCH_GetPremotor(analyfname, i, nn, bb, tbin);
            else
                decodestruct = [];
            end
            if isempty(decodestruct)
                % --- fill with nan
                decodestruct.Pdat = nan;
            end
            
            
            %% ################################# OTHER DISTANCE METRICS
            if ~isempty(timewindow)
                % ========= 1) CORR OF SMOOTHED FR
                
                % ===================================== data
                % ----------- 1) CONCATENATE ALL TRIALS (FR)
                numclasses = length(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR.classnum);
                FRmatall = [];
                Xall = [];
                IndList = cell(1,numclasses);
                count = 1;
                for cc = 1:numclasses
                    frmat = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR.classnum(cc).FRsmooth_rate_CommonTrialDur;
                    x = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR.classnum(cc).FRsmooth_xbin_CommonTrialDur;
                    if isempty(frmat)
                        continue
                    end
                    
                    if ~isempty(FRmatall) & (size(FRmatall,1) ~= size(frmat,1))
                        maxind = min([size(FRmatall,1) size(frmat,1)]);
                        FRmatall = FRmatall(1:maxind,:);
                        frmat = frmat(1:maxind,:);
                        x = x(1:maxind);
                        Xall = Xall(1:maxind,:);
                    end
                    
                    FRmatall = [FRmatall frmat];
                    Xall = [Xall x];
                    IndList{cc} = [count:count+size(frmat,2)-1];
                    count = count+size(frmat,2);
                end
                x = median(Xall,2);
                % -------------- 2) get mean fr from concatenated frmat
                FRallclass = [];
                for cc = 1:numclasses
                    frmat = FRmatall(:,IndList{cc});
                    fr = mean(frmat,2);
                    FRallclass = [FRallclass; fr'];
                end
                
                
                % -------------- 3) compute all pairwise correlations in premotor window
                %                 x = [ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR.classnum.FRsmooth_xbin_CommonTrialDur];
                %                 x = x(:,1);
                motifpredur = ALLBRANCH.alignpos(apos).ParamsFirstIter.motifpredur;
                indwind = x>=(motifpredur+timewindow(1)) & x<=(motifpredur+timewindow(2));
                RhoAll = [];
                for cc =1:numclasses
                    for ccc=cc+1:numclasses
                        
                        fr1 = FRallclass(cc,indwind);
                        fr2 = FRallclass(ccc,indwind);
                        
                        if any(isnan((fr1))) | any(isnan((fr2)))
                            continue
                        end
                        
                        rho = corr(fr1',fr2');
                        RhoAll = [RhoAll; rho];
                    end
                end
                y_rho_DAT = mean(RhoAll);
                
                
                % =========================================== Neg control
                % (shuffle)
                Nshuff = 20;
                nsamps = size(FRmatall,2);
                y_rho_shuff = [];
                for ss = 1:Nshuff
                    
                    % --- shuffle
                    indshuff = randperm(nsamps);
                    FRmatshuff = FRmatall(:, indshuff);
                    
                    % --------- 2) get mean FR
                    FRallclassSHUFF = [];
                    for cc = 1:numclasses
                        frmat = FRmatshuff(:,IndList{cc});
                        fr = mean(frmat,2);
                        FRallclassSHUFF = [FRallclassSHUFF; fr'];
                    end
                    
                    
                    % -------------- 3) compute all pairwise correlations in premotor window
                    RhoAllSHUFF = [];
                    for cc =1:numclasses
                        for ccc=cc+1:numclasses
                            
                            fr1 = FRallclassSHUFF(cc,indwind);
                            fr2 = FRallclassSHUFF(ccc,indwind);
                            
                            if any(isnan((fr1))) | any(isnan((fr2)))
                                continue
                            end
                            
                            
                            rho = corr(fr1',fr2');
                            RhoAllSHUFF = [RhoAllSHUFF; rho];
                        end
                    end
                    y_rho = mean(RhoAllSHUFF);
                    
                    y_rho_shuff = [y_rho_shuff; y_rho];
                end
                y_rho_neg = median(y_rho_shuff);
                
                % ============================================= POS CONTROL
                %             assert(numclasses == ...
                %                 length(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR_POSCONTR.classnum), 'safd');
                if numclasses ~= length(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR_POSCONTR.classnum)
                    disp('FHAIOSHFIOASHIOPDH');
                end
                
                
                % ----------- 1) CONCATENATE ALL TRIALS (FR)
                numclasses = length(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR_POSCONTR.classnum);
                FRmatall = [];
                IndList = cell(1,numclasses);
                count = 1;
                for cc = 1:numclasses
                    frmat = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR_POSCONTR.classnum(cc).FRsmooth_rate_CommonTrialDur;
                    FRmatall = [FRmatall frmat];
                    IndList{cc} = [count:count+size(frmat,2)-1];
                    count = count+size(frmat,2);
                end
                
                % -------------- 2) get mean fr from concatenated frmat
                FRallclass = [];
                for cc = 1:numclasses
                    frmat = FRmatall(:,IndList{cc});
                    fr = mean(frmat,2);
                    FRallclass = [FRallclass; fr'];
                end
                
                
                % -------------- 3) compute all pairwise correlations in premotor window
                RhoAll = [];
                for cc =1:numclasses
                    for ccc=cc+1:numclasses
                        
                        fr1 = FRallclass(cc,indwind);
                        fr2 = FRallclass(ccc,indwind);
                        
                        if any(isnan((fr1))) | any(isnan((fr2)))
                            continue
                        end
                        
                        
                        rho = corr(fr1',fr2');
                        RhoAll = [RhoAll; rho];
                    end
                end
                y_rho_POS = mean(RhoAll);
            else
                y_rho_DAT = nan;
                y_rho_neg = nan;
                y_rho_POS = nan;
            end
            
            %% get fr mat
          frinds = find(~cellfun('isempty', ...
              {ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR.classnum.FRsmooth_rate_CommonTrialDur})); % which ones are not empty
            
          frdat = struct;
          count = 1;
          for indthis = frinds
              dat = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR.classnum(indthis).FRsmooth_rate_CommonTrialDur; % extract those not empty
              frdat.classnum(count).frmat = dat;
              frdat.classnum(count).ind_classorig = indthis;
              count = count+1;
          end
%           frdat = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR.classnum(frinds); % extract those not empty
            
            
            
            %% get ff by trial
            classindsinorder = [frdat.classnum.ind_classorig];
                ffstruct = lt_neural_v2_CTXT_BRANCH_GetFF(analyfname, i, ...
                    nn, bb, prms, classindsinorder);
                
                % =========== sanity check
                % make sure num contexts match
                n1 = length(ffstruct.classnum);
                n2 = length(frdat.classnum);
                %                 n2 = size(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).ConfMatAll,1);
                %                 n2 = length(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR.classnum);
                assert(n1 ==n2, 'diff num classes ?');
                % make sure samle sizes match

                for ctocheck = 1:length(ffstruct.classnum)
                    %                     n1= size(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).FR.classnum(ctocheck).FRsmooth_rate_CommonTrialDur, 2);
                    n1=  size(frdat.classnum(ctocheck).frmat,2);
                    n2 = size(ffstruct.classnum(ctocheck).t_ff,1);
                    assert(n1 ==n2, 'diff sample size ...');
                end
            
            %% ============== save
            ind = strcmp(ListOfBranches, thisbranch);
            DatStructByBranch.bird(i).branchID(ind).regexpstr = thisbranch;
            
            DatStructByBranch.bird(i).branchID(ind).DAT.xdecode = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.xdecode; ...
                x_decode];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.ydecode= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.ydecode; ...
                y_decode];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.y_decode_z= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.y_decode_z; ...
                y_decode_z];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.ydecode_neg= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.ydecode_neg; ...
                y_decode_neg];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.ydecode_pos= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.ydecode_pos; ...
                y_decode_pos];
            
            
            DatStructByBranch.bird(i).branchID(ind).DAT.y_decode_negshuff_mean= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.y_decode_negshuff_mean; ...
                y_decode_negshuff_mean];
            DatStructByBranch.bird(i).branchID(ind).DAT.y_decode_negshuff_std= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.y_decode_negshuff_std; ...
                y_decode_negshuff_std];
            
            
            % ----------------- correlation in premotor window (decode
            % metric)
            DatStructByBranch.bird(i).branchID(ind).DAT.FdRcorr_dat = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.FRcorr_dat; ...
                y_rho_DAT];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.FRcorr_neg = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.FRcorr_neg; ...
                y_rho_neg];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.FRcorr_pos = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.FRcorr_pos; ...
                y_rho_POS];
            
            %             try
            %             DatStructByBranch.bird(i).branchID(ind).DAT.dprime= ...
            %                 [DatStructByBranch.bird(i).branchID(ind).DAT.dprime; ...
            %                 dprime];
            %             DatStructByBranch.bird(i).branchID(ind).DAT.dprime_neg= ...
            %                 [DatStructByBranch.bird(i).branchID(ind).DAT.dprime_neg; ...
            %                 dprime_neg];
            %             catch err
            %             end
            %
            
            DatStructByBranch.bird(i).branchID(ind).DAT.sylcontour_mean= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.sylcontour_mean; ...
                sylcontour_mean];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.sylcontour_x= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.sylcontour_x; ...
                sylcontour_x];
            
            location = ALLBRANCH.SummaryStruct.birds(i).neurons(nn).NOTE_Location;
            DatStructByBranch.bird(i).branchID(ind).DAT.brainregion = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.brainregion; ...
                location];
            
            % ------------ premotor decode
            DatStructByBranch.bird(i).branchID(ind).DAT.PREMOTORDECODE_pval= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.PREMOTORDECODE_pval; ...
                decodestruct.Pdat];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.PREMOTORDECODE_struct= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.PREMOTORDECODE_struct; ...
                decodestruct];
            
            
            % ---------- ff structu
            DatStructByBranch.bird(i).branchID(ind).DAT.FFstruct= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.FFstruct; ...
                ffstruct];

            DatStructByBranch.bird(i).branchID(ind).DAT.frdat= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.frdat; ...
                frdat];

            
            % ----------- inds to refer back to ALLBRANCH, SUMMARYstruct,
            % ...
            DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_apos = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_apos; ...
                apos];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_bird = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_bird; ...
                i];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_branch = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_branch; ...
                bb];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_neuron = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_neuron; ...
                nn];
            
            
            
        end
    end
end