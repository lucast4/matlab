function [AllBranch_IntDiff, AllBranch_SlopeDiff, AllBranch_SlopeOverall, ...
    AllBranch_IntCoeff_MedAbs, AllBranch_SlopeCoeff_MedAbs, AllBranch_SlopeOverallCoeff_MedAbs, ...
    AllBranch_birdnum, AllBranch_branchnum, AllBranch_neurnum] = ...
    lt_neural_v2_CTXT_Acoustic_CorrSub1(AllBirdnum, AllBranchnum, ...
    AllNeurnum, AllClassnum, AllFRmeans, AllPitch, pitch_as_predictor, ...
    SummaryStruct, doshuff)
% doshuff = 1; % then shuffles residuals, while maintaining mean pitch and
% frate for each context: i.e. null hypothesis if there is no interaction
% between context and slope.


%%
numbirds = max(AllBirdnum);
maxbranch = max(AllBranchnum);
maxneur = max(AllNeurnum);

%%

AllBranch_IntDiff = [];
AllBranch_SlopeDiff = [];
AllBranch_SlopeOverall = [];
AllBranch_birdnum = [];
AllBranch_branchnum = [];
AllBranch_neurnum = [];

AllBranch_IntCoeff_MedAbs = [];
AllBranch_SlopeCoeff_MedAbs = [];
AllBranch_SlopeOverallCoeff_MedAbs = [];

%%
for i=1:numbirds
    
    birdname = SummaryStruct.birds(i).birdname;
    
    for bb = 1:maxbranch
        
        for nn = 1:maxneur
            
            inds = AllBirdnum==i & AllBranchnum==bb & AllNeurnum==nn;
            
            if ~any(inds)
                continue
            end
            disp([num2str(i) '-' num2str(bb) '-' num2str(nn)]);
            classthis = AllClassnum{inds}';
            fratemeanthis = double(sqrt(AllFRmeans{inds}))';
            pitchthis = double(AllPitch{inds})';
            % -- center the pitch
            %             pitchthis = pitchthis - mean(pitchthis);
            % -- make class categorical
            %             classthis = nominal(classthis);
            
            if all(isnan(pitchthis))
                continue
            end
            
            % =================== do all pairs of classes
            classlist = unique(classthis)';
            CoeffAll_Intdiff = [];
            PvalAll_Intdiff = [];
            CoeffAll_Slopediff = [];
            PvalAll_Slopediff = [];
            CoeffAll_SlopeOverall = [];
            PvalAll_SlopeOverall = [];
            for c=1:length(classlist)
                for cc=c+1:length(classlist)
                    
                    class1 = classlist(c);
                    class2 = classlist(cc);
                    
                    % ------------- get data for these 2 classes only
                    indstmp = classthis==class1 | classthis==class2;
                    classtmp = classthis(indstmp);
                    fratetmp = fratemeanthis(indstmp);
                    pitchtmp = pitchthis(indstmp);
                    
                    % --- center frate and pitch
                    fratetmp = fratetmp - mean(fratetmp);
                    pitchtmp = pitchtmp - mean(pitchtmp);
                    
                    
                    %% ======================= SHUFFLE RESIDUALS?
                    if doshuff==1
                        % ====== collect all residuals (pitch, fr) pairs.
                        fratemeans = grpstats(fratetmp, classtmp);
                        pitchmeans = grpstats(pitchtmp, classtmp);
                        
                        classmeans = unique(classtmp);
                        
                        % --- make array with repeating class means
                        [~, IA] = ismember(classtmp, classmeans);
                        
                        % 1) frate
                        fratemeans_vec = fratemeans(IA);
                        frateresid_vec = fratetmp - fratemeans_vec;
                        
                        % 2) pitch
                        pitchmeans_vec = pitchmeans(IA);
                        pitchresid_vec = pitchtmp - pitchmeans_vec;
                        
                        % ------------ shuffle trials, while keeping resid
                        % for pitch and frate matched
                        indshuff = randperm(length(fratetmp));
                        
                        fratetmp_SHUFF = fratemeans_vec + frateresid_vec(indshuff);
                        pitchtmp_SHUFF = pitchmeans_vec + pitchresid_vec(indshuff);
                        
                        % ================= reassign back to data
                        fratetmp = fratetmp_SHUFF;
                        pitchtmp = pitchtmp_SHUFF;
                        
                        if (0)
                            % === plot, sanity check
                            lt_figure; hold on;
                            % --- dat
                            lt_subplot(2,2,1); hold on;
                            title('dat');
                            plot(pitchtmp(classtmp==class1), fratetmp(classtmp==class1), 'ob');
                            plot(pitchtmp(classtmp==class2), fratetmp(classtmp==class2), 'or');
                            % --- shuff
                            lt_subplot(2,2,2); hold on;
                            title('SHUFF');
                            plot(pitchtmp_SHUFF(classtmp==class1), ...
                                fratetmp_SHUFF(classtmp==class1), 'ob');
                            plot(pitchtmp_SHUFF(classtmp==class2), ...
                                fratetmp_SHUFF(classtmp==class2), 'or');
                        end
                        
                    end
                    
                    %% continue - do regression
                    % --- make class categorical
                    classtmp = categorical(classtmp);
                    
                    % --- zscore everything
                    fratetmp = (fratetmp - mean(fratetmp))./std(fratetmp);
                    pitchtmp = (pitchtmp - mean(pitchtmp))./std(pitchtmp);
                    
                    % ----------- fit linear model
                    tbl = table(classtmp, fratetmp, pitchtmp);
                    %             tbl = table(classthis(indstmp), fratemeanthis(indstmp), pitchthis(indstmp), 'VariableNames', ...
                    %                 {'classthis', 'fratemeanthis', 'pitchthis'});
                    %             tbl = table(classthis, fratemeanthis, pitchthis, 'VariableNames', ...
                    %                 {'classthis', 'fratemeanthis', 'pitchthis'});
                    if pitch_as_predictor==1
                        %  then pitch as predictor
                        formula = 'fratetmp ~ pitchtmp*classtmp';
                    else
                        %  then frate as predictor
                        formula = 'pitchtmp ~ fratetmp*classtmp';
                    end
                    lme = fitlme(tbl, formula);
                    
                    % =========== collect effects and pvals
                    % --- intercept difference
                    tmp = string(classtmp(end));
                    strtofind = ['classtmp_' tmp{1}];
                    coeff_intdiff = lme.Coefficients.Estimate(strcmp(lme.CoefficientNames, strtofind));
                    pval_intdiff = lme.Coefficients.pValue(strcmp(lme.CoefficientNames, strtofind));
                    assert(length(coeff_intdiff)==1, 'asdfas');
                    
                    % --- slope difference
                    if pitch_as_predictor==1
                        strtofind2 = [strtofind ':pitchtmp'];
                    else
                        strtofind2 = [strtofind ':fratetmp'];
                    end
                    coeff_slopediff = lme.Coefficients.Estimate(strcmp(lme.CoefficientNames, strtofind2));
                    pval_slopediff = lme.Coefficients.pValue(strcmp(lme.CoefficientNames, strtofind2));
                    assert(length(coeff_slopediff)==1,'asdfas');
                    
                    % ---- signficant overall slope
                    if pitch_as_predictor==1
                        strtofind3 = 'pitchtmp';
                    else
                        strtofind3 = 'fratetmp';
                    end
                    coeff_slopeOverall = lme.Coefficients.Estimate(strcmp(lme.CoefficientNames, strtofind3));
                    pval_slopeOverall = lme.Coefficients.pValue(strcmp(lme.CoefficientNames, strtofind3));
                    assert(length(coeff_slopeOverall)==1,'asdfas');
                    
                    
                    % ------------------------ SAVE
                    CoeffAll_Intdiff = [CoeffAll_Intdiff; coeff_intdiff];
                    PvalAll_Intdiff = [PvalAll_Intdiff; pval_intdiff];
                    CoeffAll_Slopediff = [CoeffAll_Slopediff; coeff_slopediff];
                    PvalAll_Slopediff = [PvalAll_Slopediff; pval_slopediff];
                    
                    CoeffAll_SlopeOverall = [CoeffAll_SlopeOverall; coeff_slopeOverall];
                    PvalAll_SlopeOverall = [PvalAll_SlopeOverall; pval_slopeOverall];
                end
            end
            
            % ================ ARE ANY SLOPE OR INTERCEPT DIFFERENT BETWEEN
            % ANY PAIR?
            % ----- Do Bonferonni correction on pvals
            PvalAll_Intdiff = PvalAll_Intdiff.*length(PvalAll_Intdiff);
            PvalAll_Slopediff = PvalAll_Slopediff.*length(PvalAll_Slopediff);
            PvalAll_SlopeOverall = PvalAll_SlopeOverall.*length(PvalAll_SlopeOverall);
            
            inteffect = any(PvalAll_Intdiff<0.05);
            slopeeffect = any(PvalAll_Slopediff<0.05);
            slopeoveralleffect = any(PvalAll_SlopeOverall<0.05);
            
            % ------------------ COLLECT MEDIANS OF ABSOLUTE VALUES OF THE COEFFICIENTS
            intcoeff = median(abs(CoeffAll_Intdiff));
            slopecoeff = median(abs(CoeffAll_Slopediff));
            slopeoverallcoeff = median(abs(CoeffAll_SlopeOverall));
            
            % =================== OUTPUT
            AllBranch_IntDiff = [AllBranch_IntDiff; inteffect];
            AllBranch_SlopeDiff = [AllBranch_SlopeDiff; slopeeffect];
            AllBranch_SlopeOverall = [AllBranch_SlopeOverall; slopeoveralleffect];
            
            AllBranch_IntCoeff_MedAbs = [AllBranch_IntCoeff_MedAbs; intcoeff];
            AllBranch_SlopeCoeff_MedAbs = [AllBranch_SlopeCoeff_MedAbs; slopecoeff];
            AllBranch_SlopeOverallCoeff_MedAbs = [AllBranch_SlopeOverallCoeff_MedAbs; slopeoverallcoeff];
            
            AllBranch_birdnum = [AllBranch_birdnum; i];
            AllBranch_branchnum = [AllBranch_branchnum; bb];
            AllBranch_neurnum = [AllBranch_neurnum; nn];
            
            if (0)
                %%
                lt_figure; hold on;
                %                 title([birdname '-' branchname{1} '-n' num2str(nn)]);
                ylabel('pitch');
                xlabel('mean frate (sqrt)');
                assert(sum(inds)==1,'asdfsd');
                
                classthis = AllClassnum{inds}';
                fratemeanthis = double(sqrt(AllFRmeans{inds}))';
                pitchthis = double(AllPitch{inds})';
                
                % ---------------------------
                pcols = lt_make_plot_colors(max(classthis), 0,0);
                for cc=1:max(classthis)
                    indtmp = classthis==cc;
                    
                    x = fratemeanthis(indtmp);
                    y = pitchthis(indtmp);
                    
                    plot(x,y, 'x', 'Color', pcols{cc});
                    
                    % -------- plot mean pitch and frate for this class
                    xmean = mean(x);
                    xsem = lt_sem(x);
                    ymean = mean(y);
                    ysem = lt_sem(y);
                    
                    %                lt_plot(xmean, ymean, {'Errors', ysem, 'Color', pcols{cc}});
                    lt_plot(xmean, ymean, {'Color', pcols{cc}});
                    line(xmean+[-xsem xsem], [ymean ymean], 'Color', pcols{cc});
                    line([xmean xmean], ymean+[-ysem ysem], 'Color', pcols{cc});
                end
                
            end
        end
    end
end
