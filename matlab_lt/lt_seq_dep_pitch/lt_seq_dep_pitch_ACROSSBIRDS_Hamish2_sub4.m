function [OUT_Generalization_Nontarg, OUT_CVchange_TargNontarg] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub4(DATSTRUCT)
%% lt 7/27/18 - ALL PLOTS LOOKING AT TARGET SYL + same type ONLY

%%

% ====== FIRST, subsample all to only iclude target syls
inds = find(DATSTRUCT.All_Issame==1);
DATSTRUCT = lt_structure_subsample_all_fields(DATSTRUCT, inds, 1);


% ===== SECOND, extract all fields to qworkspace variables
fnames = fieldnames(DATSTRUCT)';
for fn=fnames
    fn = fn{1};
    eval([fn ' = DATSTRUCT.(fn);'])
end

%% ########## RELATIONSHIP BETWEEN AFP BIAS AT TARG AND GENERALIZATION

% ====================== 1) generalization as a function of baseline AFP
% bias
lt_figure; hold on;

% ========================= 1) plot learning for all targets
lt_subplot(3,2,1); hold on;
title('learning');
ylabel('learning at targ (PBS minus PBS)');
xlabel('baseline AFP bias at target (in direction of laerning)');

ind = All_Istarg==1; assert(sum(ind)== max(All_exptcounter))
x = All_learndir(ind).*(All_FF_BASE_PBS(ind) - All_FF_BASE_MUSC(ind));
y_targ = All_learndir(ind).*(All_FF_WN_PBS(ind) - All_FF_BASE_PBS(ind));

plot(x, y_targ, 'ok');
lt_plot_zeroline;


% ==========================
lt_subplot(3,2,2); hold on;
title('generalization');
ylabel('generalization (PBS minus PBS)');
xlabel('baseline AFP bias at target (in direction of laerning)');

numexpts = max(All_exptcounter);

x_targ_all = nan(numexpts,1);
y_targ_all = nan(numexpts,1);
cvdiff_PBS = nan(numexpts,1);
cvdiff_MUSC = nan(numexpts,1);

y_nontarg_all = cell(numexpts,1);
cvdiff_PBS_nontarg = cell(numexpts,1);
cvdiff_MUSC_nontarg = cell(numexpts,1);
ffdiff_MUSC_nontarg = cell(numexpts,1);
ffdiff_PBS_nontarg = cell(numexpts,1);

for i=1:numexpts
    
    % -- targ
    ind = All_exptcounter==i & All_Istarg==1; assert(sum(ind)==1);
    x = All_learndir(ind).*(All_FF_BASE_PBS(ind) - All_FF_BASE_MUSC(ind));
    y_targ = All_learndir(ind).*(All_FF_WN_PBS(ind) - All_FF_BASE_PBS(ind));
    
    x_targ_all(i) = x;
    y_targ_all(i) = y_targ;
    cvdiff_PBS(i) = All_CV_WN_PBS(ind) - All_CV_BASE_PBS(ind);
    cvdiff_MUSC(i) = All_CV_WN_MUSC(ind) - All_CV_BASE_MUSC(ind);
    
    % -- nontarg
    ind = All_exptcounter==i & All_Istarg==0; assert(length(unique(All_learndir(ind)))<2);
    Y = All_learndir(ind).*(All_FF_WN_PBS(ind) - All_FF_BASE_PBS(ind));
    Y = Y./y_targ;
    y_nontarg_all{i} = Y;
    
    cvdiff_PBS_nontarg{i} = All_CV_WN_PBS(ind) - All_CV_BASE_PBS(ind);
    cvdiff_MUSC_nontarg{i} = All_CV_WN_MUSC(ind) - All_CV_BASE_MUSC(ind);
    
    ffdiff_MUSC_nontarg{i} = All_learndir(ind).*(All_FF_WN_MUSC(ind) - All_FF_BASE_MUSC(ind));
    ffdiff_PBS_nontarg{i} = All_learndir(ind).*(All_FF_WN_PBS(ind) - All_FF_BASE_PBS(ind));
    
    if ~isempty(Y)
        plot(x, Y, 'ok');
    end
    
    if length(Y)>1
        line([x x], [min(Y) max(Y)], 'Color', [0.6 0.6 0.6]);
    end
end

OUT_Generalization_Nontarg = y_nontarg_all;

% =============
lt_subplot(3,2,3); hold on;
xlabel('baseline AFP bias at target (in direction of laerning)');
ylabel('generalization (PBS minus PBS)');

x = x_targ_all;
y = cellfun(@mean, y_nontarg_all);
lt_regress(y, x, 1);
lt_plot_zeroline;


% =============== compare cases with bias positive vs. negative.
lt_subplot(3,2,4); hold on;
xlabel('AGAINST -- TOWARDS');
ylabel('generalization (fraction)');

Y = {};
% -- against
indtmp = x_targ_all<0;
ythis = cellfun(@mean, y_nontarg_all(indtmp));
% ythis = cell2mat(y_nontarg_all(indtmp));
Y{1} = ythis;
% -- towards
indtmp = x_targ_all>0;
ythis = cellfun(@mean, y_nontarg_all(indtmp));
% ythis = cell2mat(y_nontarg_all(indtmp));
Y{2} = ythis;
% --- plot
lt_plot_bar([1 2]+0.15, cellfun(@mean, Y), {'Errors', cellfun(@lt_sem, Y), 'Color', 'r'});
lt_plot_MultDist(Y, [1 2], 0, 'k', 1);
lt_plot_zeroline;
[~, p] = ttest2(Y{1}, Y{2});
lt_plot_pvalue(p, 'ttest', 1);

OUT_Generalization_AgainstTowards = Y;

% ===============
lt_subplot(3,2,5); hold on;
xlabel('AGAINST -- TOWARDS');
title('laerning at targ (hz)');

Y = {};
% -- against
indtmp = x_targ_all<0;
ythis = y_targ_all(indtmp);
Y{1} = ythis;
% -- towards
indtmp = x_targ_all>0;
ythis = y_targ_all(indtmp);
Y{2} = ythis;
% --- plot
lt_plot_bar([1 2]+0.15, cellfun(@mean, Y), {'Errors', cellfun(@lt_sem, Y), 'Color', 'r'});
lt_plot_MultDist(Y, [1 2], 0, 'k', 1);
lt_plot_zeroline;
[~, p] = ttest2(Y{1}, Y{2});
lt_plot_pvalue(p, 'ttest', 1);


% ===================
lt_subplot(3,2,6); hold on;
xlabel('targ afp bias (dir of laerning)');
ylabel('nontarg learn (bk = PBS, rd=MUSC)');

for i=1:length(x_targ_all)
    x = x_targ_all(i);
    
    yPBS = mean(ffdiff_PBS_nontarg{i});
    yMUSC = mean(ffdiff_MUSC_nontarg{i});
    
    plot(x, yPBS, 'ok');
    plot(x, yMUSC, 'or');
    %    line([x x], [yPBS yMUSC], 'Color', [0.6 0.6 0.6]);
end

%  TO DO:
% =============== generalization as function of 1) amount of learning, 2)
% amount of learning consolidated




%% ====== CV
lt_figure; hold on;


% ============= 1) CV for targ
lt_subplot(3,2,1); hold on;
xlabel('base AFP bias, target (dir of learning)');
ylabel('change in CV');
title('[TARGET] PBS=k; MUSC=r');
% ---- PBS
x = x_targ_all;
y = cvdiff_PBS;
plot(x,y, 'ok');

% ---- MUSC
x = x_targ_all;
y = cvdiff_MUSC;
plot(x,y, 'or');

% --
lt_plot_zeroline;


% ============ 2) CV of nontarg [all nontarg]
lt_subplot(3,2,2); hold on;
xlabel('base AFP bias, target (dir of learning)');
ylabel('change in CV (NONTARG)');
title('[NONTARG] PBS=k; MUSC=r');

% --- PBS
pcol = 'k';
for i=1:length(x_targ_all)
    plot(x_targ_all(i), cvdiff_PBS_nontarg{i}, 'o', 'Color', pcol);
end

% --- MUSC
pcol = 'r';
for i=1:length(x_targ_all)
    plot(x_targ_all(i), cvdiff_MUSC_nontarg{i}, 'o', 'Color', pcol);
end

lt_plot_zeroline;


% ============ 2) CV of nontarg [mean across nontargs];
lt_subplot(3,2,3); hold on;
xlabel('base AFP bias, target (dir of learning)');
ylabel('change in CV (NONTARG)');
title('[NONTARG] PBS=k; MUSC=r');

% --- PBS
pcol = 'k';
x = x_targ_all;
y = cellfun(@mean, cvdiff_PBS_nontarg);
plot(x,y, 'o', 'Color', pcol);

% --- MUSC
pcol = 'r';
x = x_targ_all;
y = cellfun(@mean, cvdiff_MUSC_nontarg);
plot(x,y, 'o', 'Color', pcol);

lt_plot_zeroline;


% ============ 3) BAR PLOTS
lt_subplot(3,2,4); hold on;
title('PBS');
xlabel('AGAINST(targ) AGAINST (same) -- TOWARDS(targ) TOWARDS(same)');
ylabel('change in CV');

% ------------ PBS
Y = [cvdiff_PBS cellfun(@mean, cvdiff_PBS_nontarg)];

% against
indstmp = x_targ_all<0;
xthis = [1 2];
Ythis = Y(indstmp,:);

plot(xthis, Ythis, '-k');
ymean = mean(Ythis);
ysem = lt_sem(Ythis);
lt_plot_bar(xthis+0.15, ymean, {'Errors', ysem});
[~, p] = ttest(Ythis(:,1), Ythis(:,2));
if p<0.1
    lt_plot_text(mean(xthis), 0.015, ['p=' num2str(p)], 'r');
end

% towards
indstmp = x_targ_all>0;
xthis = [5 6];
Ythis = Y(indstmp,:);

plot(xthis, Ythis, '-k');
ymean = mean(Ythis);
ysem = lt_sem(Ythis);
lt_plot_bar(xthis+0.15, ymean, {'Errors', ysem});
[~, p] = ttest(Ythis(:,1), Ythis(:,2));
if p<0.1
    lt_plot_text(mean(xthis), 0.015, ['p=' num2str(p)], 'r');
end

lt_plot_zeroline;
OUT_CVchange_TargNontarg = Y;

% ============ 3) BAR PLOTS
lt_subplot(3,2,5); hold on;
title('MUSC');
xlabel('AGAINST(targ) AGAINST (same) -- TOWARDS(targ) TOWARDS(same)');
ylabel('change in CV');

% ------------ PBS
Y = [cvdiff_MUSC cellfun(@mean, cvdiff_MUSC_nontarg)];

% against
indstmp = x_targ_all<0;
xthis = [1 2];
Ythis = Y(indstmp,:);

plot(xthis, Ythis, '-k');
ymean = mean(Ythis);
ysem = lt_sem(Ythis);
lt_plot_bar(xthis+0.15, ymean, {'Errors', ysem});
[~, p] = ttest(Ythis(:,1), Ythis(:,2));
if p<0.1
    lt_plot_text(mean(xthis), 0.015, ['p=' num2str(p)], 'r');
end

% towards
indstmp = x_targ_all>0;
xthis = [5 6];
Ythis = Y(indstmp,:);

plot(xthis, Ythis, '-k');
ymean = mean(Ythis);
ysem = lt_sem(Ythis);
lt_plot_bar(xthis+0.15, ymean, {'Errors', ysem});
[~, p] = ttest(Ythis(:,1), Ythis(:,2));
if p<0.1
    lt_plot_text(mean(xthis), 0.015, ['p=' num2str(p)], 'r');
end

lt_plot_zeroline;


%% ######################## USING BASELINE BIAS AT NONTARG SYLS
lt_figure; hold on

% ==========================
lt_subplot(3,2,1); hold on;
xlabel('nontarg AFP bias (baseline, targ learn dir)');
ylabel('CV change (PBS)');

indstmp = All_Istarg==0;
x = All_learndir(indstmp).*(All_FF_BASE_PBS(indstmp) - All_FF_BASE_MUSC(indstmp));
y = All_CV_WN_PBS(indstmp) - All_CV_BASE_PBS(indstmp);

plot(x,y,'ok');
lt_plot_zeroline;

% ====================
lt_subplot(3,2,2); hold on;
xlabel('nontarget bias vs. learning: AGAINST -- TOWARDS');
ylabel('CV change (PBS)');

Y = {};
X = [1 2];
% --- against
Y{1} = y(x<0);
% --- towards
Y{2} = y(x>0);
% ------------ PLOT
lt_plot_MultDist(Y, X, 0, 'k', 1);
lt_plot_bar(X, cellfun(@mean, Y), {'Errors', cellfun(@lt_sem, Y), 'Color', 'r'});


% ==========================
lt_subplot(3,2,3); hold on;
xlabel('nontarg AFP bias (baseline, targ learn dir)');
ylabel('learning (PBS, hz)');

indstmp = All_Istarg==0;
x = All_learndir(indstmp).*(All_FF_BASE_PBS(indstmp) - All_FF_BASE_MUSC(indstmp));
y =  All_learndir(indstmp).*(All_FF_WN_PBS(indstmp) - All_FF_BASE_PBS(indstmp));

plot(x,y,'ok');
lt_plot_zeroline;


% ====================
lt_subplot(3,2,4); hold on;
xlabel('nontarget bias vs. learning: AGAINST -- TOWARDS');
ylabel('Learn (PBS, hz)');

Y = {};
X = [1 2];
% --- against
Y{1} = y(x<0);
% --- towards
Y{2} = y(x>0);
% ------------ PLOT
lt_plot_MultDist(Y, X, 0, 'k', 1);
lt_plot_bar(X, cellfun(@mean, Y), {'Errors', cellfun(@lt_sem, Y), 'Color', 'r'});








