function lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub6(DATSTRUCT)

%% PLOTS pitch contour related stuff


%% extract variables

fnames = fieldnames(DATSTRUCT)';
for fn=fnames
    fn = fn{1};
    eval([fn ' = DATSTRUCT.(fn);'])
end


%% ####################### PLOT SUMMARY [wiggles, across experiments]
lt_figure; hold on;

% ========= HI minus LO
lt_subplot(3,2,1); hold on;
xlabel('base bias');
ylabel('wiggle (HI) - wiggle (LO)');
title('wiggle = median of std');

indsgood = ~isnan(BaseBiasAll); % i.e. smalle sample size some

% wigglemat = cellfun(@median, WiggleAll, 'UniformOutput', 0);
wigglemat = cellfun(@mean, WiggleAll, 'UniformOutput', 0);
wigglemat = cell2mat(wigglemat(indsgood,:));
biasmat = BaseBiasAll(indsgood);

plot(biasmat, wigglemat(:,2) - wigglemat(:,1), 'ok')
lt_regress(wigglemat(:,2) - wigglemat(:,1), biasmat, 1);
lt_plot_zeroline;


% ============== wiggle vs. ff, one corr for each syl
lt_subplot(3,2,2); hold on;
xlabel('base AFP bias');
ylabel('trial corr (ffmean vs. wiggle)');

indsgood = ~isnan(BaseBiasAll); % i.e. smalle sample size some
x = BaseBiasAll(indsgood);
y = WiggleVsMeanFF(indsgood);
lt_regress(y, x, 1);
lt_plot_zeroline;


% ================= SEPARATE BY DIRECTION OF BASELINE AFP BIAS
lt_subplot(3,2,4); hold on;

Y = {};

% -- neg bias
indstmp = BaseBiasAll<0;
wigglediff = cell2mat(WiggleAll(indstmp,2)) - cell2mat(WiggleAll(indstmp,1));
Y{1} = wigglediff;

% -- pos bias
indstmp = BaseBiasAll>0;
wigglediff = cell2mat(WiggleAll(indstmp,2)) - cell2mat(WiggleAll(indstmp,1));
Y{2} = wigglediff;

% -- plot
lt_plot_MultDist(Y, [1 2], 1, 'k', 1);
xlim([0 3]);
xlabel('BaseAFPbias (DN -- UP)');
ylabel('wiggle(HIrends) - wiggle(LO)');



% ========= HI minus LO
lt_subplot(3,2,5); hold on;
xlabel('base bias');
ylabel('wiggle (HI) - wiggle (LO)');
title('ONLY IF WIGGLE > MUSC');

indsgood = ~isnan(BaseBiasAll) & DATSTRUCT.WiggleOverMusc==1; % i.e. smalle sample size some


% wigglemat = cellfun(@median, WiggleAll, 'UniformOutput', 0);
wigglemat = cellfun(@mean, WiggleAll, 'UniformOutput', 0);
wigglemat = cell2mat(wigglemat(indsgood,:));
biasmat = BaseBiasAll(indsgood);

plot(biasmat, wigglemat(:,2) - wigglemat(:,1), 'ok')
lt_regress(wigglemat(:,2) - wigglemat(:,1), biasmat, 1);
lt_plot_zeroline;


% ======== PLOT FRACTION OF CASES THAT HAVE GREATER WIGGLE THAN MUSC
lt_subplot(3, 2,6); hold on;
title('wiggles (bu: POS BIAS ... rd: NEG BIAS)');
xlabel('low rends');
ylabel('high rends');

indsgood = ~cellfun(@isempty, DATSTRUCT.WiggleAll(:,1)); % not empty

wigmat = cell2mat(DATSTRUCT.WiggleAll(indsgood, :));
biasmat = DATSTRUCT.BaseBiasAll(indsgood);

% --- POSITIVE BASE BIAS
indtmp = biasmat>0;
pcol = 'b';

plot(wigmat(indtmp,1), wigmat(indtmp, 2), 'x', 'Color', pcol);
lt_plot(mean(wigmat(indtmp,1)), mean(wigmat(indtmp, 2)), {'Color', pcol});

% --- NEGATIVE BASE BIAS
indtmp = biasmat<0;
pcol = 'r';

plot(wigmat(indtmp,1), wigmat(indtmp, 2), 'x', 'Color', pcol);
lt_plot(mean(wigmat(indtmp,1)), mean(wigmat(indtmp, 2)), {'Color', pcol});

% ----
lt_plot_makesquare_plot45line(gca, 'k');




