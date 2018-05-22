%% ====== LINEAR MODEL TESTING WHETHER SLOPE DEPENDS ON BREGION

y1 = AllPairs_Means(:,1);
y2 = AllPairs_Means(:,2);
tbl = table(y1, y2, AllPairs_Birdnum, AllPairs_Bregions, 'VariableNames', ...
    {'y1', 'y2', 'birdnum', 'bregion'});

% mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum) + (-1+bregion:y1|birdnum)';  % full model
mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum)';  % full model
% mdl = 'y2 ~ y1 + bregion + bregion:y1 + (-1 + y1|birdnum) + (-1 + bregion|birdnum)';  % variable intercept by brain region
lme = fitlme(tbl, mdl)


%% =========== SAME, BUT RESTRICTED TO BIRDS WITH BOTH DATA
goodbirds = intersect(AllPairs_Birdnum(strcmp(AllPairs_Bregions, 'LMAN')) ,...
    AllPairs_Birdnum(strcmp(AllPairs_Bregions, 'RA')));
indsgood = ismember(AllPairs_Birdnum, goodbirds');

y1 = AllPairs_Means(:,1);
y2 = AllPairs_Means(:,2);
tbl = table(y1(indsgood), y2(indsgood), AllPairs_Birdnum(indsgood), AllPairs_Bregions(indsgood), 'VariableNames', ...
    {'y1', 'y2', 'birdnum', 'bregion'});

mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum) + (-1+bregion:y1|birdnum)';  % full model
lme = fitlme(tbl, mdl)

%% ============ SAME, BUT RESTRICTED TO BIRDS WITH ONLY ONE DATASET
goodbirds = intersect(AllPairs_Birdnum(strcmp(AllPairs_Bregions, 'LMAN')) ,...
    AllPairs_Birdnum(strcmp(AllPairs_Bregions, 'RA')));
indsgood = ~ismember(AllPairs_Birdnum, goodbirds');

y1 = AllPairs_Means(:,1);
y2 = AllPairs_Means(:,2);
tbl = table(y1(indsgood), y2(indsgood), AllPairs_Birdnum(indsgood), AllPairs_Bregions(indsgood), 'VariableNames', ...
    {'y1', 'y2', 'birdnum', 'bregion'});

% mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum) + (-1+bregion:y1|birdnum)';  % full model
mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum)';  % full model
lme = fitlme(tbl, mdl)

%% =========== SAME, BUT RESTRICTED TO BIRDS THAT I RECORDED
% -------- FIND THE GOOD BIRDS
goodbirds = [];
for j=1:length(SummaryStruct.birds)
    if isfield(SummaryStruct.birds(j).neurons(1), 'isRAsobermel')
        continue
    else
        % -- then is my bird
        goodbirds = [goodbirds j];
    end
end
disp({SummaryStruct.birds(goodbirds).birdname});

% ------- INDS FOR THESE BIRDS
indsgood = ismember(AllPairs_Birdnum, goodbirds);

y1 = AllPairs_Means(:,1);
y2 = AllPairs_Means(:,2);
tbl = table(y1(indsgood), y2(indsgood), AllPairs_Birdnum(indsgood), AllPairs_Bregions(indsgood), 'VariableNames', ...
    {'y1', 'y2', 'birdnum', 'bregion'});

mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum) + (-1+bregion:y1|birdnum)';  % full model
mdl = 'y2 ~ y1 + bregion:y1 + (-1 + y1|birdnum)';  % full model
lme = fitlme(tbl, mdl)

%% =============== SAME, BUT DO MODEL FOR SINGLE BIRD
birdtodo = 'wh44wh39';
goodbirds = [];
for j=1:length(SummaryStruct.birds)
    if ~strcmp(SummaryStruct.birds(j).birdname, birdtodo)
        continue
    else
        % -- then is my bird
        goodbirds = [goodbirds j];
    end
end
goodbirds = find(strcmp({SummaryStruct.birds.birdname}, birdtodo));

% ------- INDS FOR THESE BIRDS
indsgood = ismember(AllPairs_Birdnum, goodbirds);

% -----
y1 = AllPairs_Means(:,1);
y2 = AllPairs_Means(:,2);
tbl = table(y1(indsgood), y2(indsgood), AllPairs_Birdnum(indsgood), AllPairs_Bregions(indsgood), 'VariableNames', ...
    {'y1', 'y2', 'birdnum', 'bregion'});

mdl = 'y2 ~ y1 + bregion:y1';  % full model
lme = fitlme(tbl, mdl)

%% ###############################################
%% ###############################################

% ========= [PLOT] fixed effect coefficients

fixednames = lme.CoefficientNames;
fixedeffCI = lme.coefCI;
fixedeff = lme.fixedEffects;
fixedeffSE = lme.Coefficients.SE;

lt_figure; hold on;
x = 1:length(fixednames);
errorbar(x, fixedeff,  fixedeff-fixedeffCI(:,1), -fixedeff+fixedeffCI(:,2), 'LineStyle', 'none');
plot(x, fixedeff, 'ok');
xlim([0 x(end)+1]);
lt_plot_zeroline;
set(gca, 'XTick', x, 'XTickLabel', fixednames);
rotateXLabels(gca, 45);

% --------- plot p values
YLIM = ylim;
pvals = lme.Coefficients.pValue;
for x=1:length(pvals)
    if pvals(x)<0.1
   lt_plot_text(x-0.3, YLIM(2), ['p=' num2str(pvals(x))], 'r', 10)
    end
end

ylim([-1 1]);

% ========== [PLOT] BUT REPRESENTING AS LMAN/RA
lt_figure; hold on;
coeff_LMAN = fixedeff(2);
coeffCI_LMAN = fixedeffCI(2,:);
coeffSE_LMAN = fixedeffSE(2);

coeff_RA = fixedeff(2) + fixedeff(3);
coeffCI_RA = fixedeff(2) + fixedeffCI(3,:);
coeffSE_RA = fixedeffSE(3);

X = [1 2];
lt_plot_bar(X, [coeff_LMAN coeff_RA], {'Errors', [coeffSE_LMAN coeffSE_RA], 'Color', 'k'});

lt_plot_pvalue(lme.Coefficients.pValue(3), 'RA vs LMAN', 1);
%% ========= residual plot
lt_figure; hold on;
plotResiduals(lme, 'fitted');

%% =============== PLOT RANDOM EFFECTS
[B, Bnames, Stats] = randomEffects(lme);

% ========== plot all random effects
lt_figure; hold on;
ylabel('random effects');
Xnames = {};
for x=1:length(Bnames.Name)
    y = B(x);
    yCI = [Stats.Lower(x) Stats.Upper(x)];
    p = Stats.pValue(x);
    
    errorbar(x, y,  y-yCI(1), -y+yCI(2), 'LineStyle', 'none', 'Color', 'k');
    plot(x, y, 'ok');
    
    % --- pvale
    if p<0.1
        lt_plot_text(x, y, ['p=' num2str(p)], 'r', 8);
    end
    % --- collect xlabels
    Xnames = [Xnames [Bnames.Name{x} '-' Bnames.Group{x} '-' Bnames.Level{x}]];
end
set(gca, 'XTick', 1:length(Xnames), 'XTickLabel', Xnames);

lt_plot_zeroline;

rotateXLabels(gca, 45);


ylim([-1 1]);

%% ================ PLOT COEFFICIENTS FOR EACH BIRD (SUM OF FIXED AND RANDOM)
%  IN PROGRESS

coeffname = 'y1';

coeff_fixed = fixedeff(strcmp(fixednames, coeffname));

% ===== for 
intersect(fixednames, Stats.Name)


%% ================ FOR EACH BIRD OVERLAY FIT WITH DATA

% --- each bird plot data along with regression lines (and simulated lines)
x = lme.Variables;
F = fitted(lme);

figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


numbirds = max(x.birdnum);

for i=1:numbirds
  [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
  title(SummaryStruct.birds(i).birdname);
  ylabel('lines: fit');
  
  % ####################### ACTUAL DAT
  % ===== LMAN
  pcol = 'b';
  inds = x.birdnum==i & strcmp(x.bregion, 'LMAN');
  if any(inds)
  xdat = x.y1(inds);
  ydat = x.y2(inds);
  
  plot(xdat, ydat, 'o', 'Color', pcol)
  end
  
  % ==== RA
  pcol = 'r';
  inds = x.birdnum==i & strcmp(x.bregion, 'RA');
  if any(inds)
  xdat = x.y1(inds);
  ydat = x.y2(inds);
  
  plot(xdat, ydat,  'o', 'Color', pcol)
  end
  
  % ########################## FIT DAT
   % ===== LMAN
  pcol = 'b';
  inds = x.birdnum==i & strcmp(x.bregion, 'LMAN');
  if any(inds)
  xdat = x.y1(inds);
  ydat = F(inds);
  
  % ---------- plot a line
  [~, indtmp1] = min(xdat);
  [~, indtmp2] = max(xdat);
  xthis = [xdat(indtmp1) xdat(indtmp2)];
  ythis = [ydat(indtmp1) ydat(indtmp2)];
  line(xthis, ythis, 'Color', pcol, 'LineWidth', 2);
  end
  
  % ==== RA
  pcol = 'r';
  inds = x.birdnum==i & strcmp(x.bregion, 'RA');
  if any(inds)
  xdat = x.y1(inds);
  ydat = F(inds);
  
  % ---------- plot a line
  [~, indtmp1] = min(xdat);
  [~, indtmp2] = max(xdat);
  xthis = [xdat(indtmp1) xdat(indtmp2)];
  ythis = [ydat(indtmp1) ydat(indtmp2)];
  line(xthis, ythis, 'Color', pcol, 'LineWidth', 2);
  end
  
  % ==========
  axis tight;
    lt_plot_makesquare_plot45line(gca, 'k', -1);

 
end



%% 

%% ================ FOR EACH BIRD OVERLAY [PREDICTION] WITH DATA

% --- each bird plot data along with regression lines (and simulated lines)
x = lme.Variables;
y = lme.response;
F = random(lme, x);

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


numbirds = max(x.birdnum);

for i=1:numbirds
  [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
  title(SummaryStruct.birds(i).birdname);
  ylabel('solid: fit');
  
  % ####################### ACTUAL DAT
  % ===== LMAN
  pcol = 'b';
  inds = x.birdnum==i & strcmp(x.bregion, 'LMAN');
  if any(inds)
  xdat = x.y1(inds);
  ydat = x.y2(inds);
  
  plot(xdat, ydat, 'o', 'Color', pcol)
  end
  
  % ==== RA
  pcol = 'r';
  inds = x.birdnum==i & strcmp(x.bregion, 'RA');
  if any(inds)
  xdat = x.y1(inds);
  ydat = x.y2(inds);
  
  plot(xdat, ydat,  'o', 'Color', pcol)
  end
  
  % ########################## FIT DAT
   % ===== LMAN
  pcol = 'b';
  inds = x.birdnum==i & strcmp(x.bregion, 'LMAN');
  if any(inds)
  xdat = x.y1(inds);
  ydat = F(inds);
  
  lt_plot(xdat, ydat, {'Color', pcol});
  end
  
  % ==== RA
  pcol = 'r';
  inds = x.birdnum==i & strcmp(x.bregion, 'RA');
  if any(inds)
  xdat = x.y1(inds);
  ydat = F(inds);
  
  lt_plot(xdat, ydat, {'Color', pcol});
  end
  
  % ==========
    lt_plot_makesquare_plot45line(gca, 'k', -1);

 
end


%% ===== [plot] FOR EACH BIRD PLOT SLOPE (LMAN/RA)



