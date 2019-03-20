function [rhodat, rhoperm_all] = lt_neural_POPLEARN_Xcov_EpochScal_sub(allDat_xcov, allDat_learn, ...
    allbnum, allenum, allswnum, doshift, bintobinchanges, scalwind, syltype, ...
    nshuff, vers)
if ~exist('vers', 'var')
    vers = 1; % original, lookign at average xcov change
end
    
%%
if vers==1
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
elseif vers==2
   % looking not at aveage chagne in xcov, but at change in neural-ff corr (i.e. FF split, AFp bias..)
   X = squeeze(allDat_xcov);
   X = X(2:end,:); % remove row of 0, since wil add in subsequent code.
else
    sdfasdf
end
Y = squeeze(allDat_learn(:,:, syltype, :));
[exptID, exptID_U] = lt_tools_grp2idx({allbnum, allenum, allswnum});

% --- add 0, since they are all deviation from baseline
Xneural = [zeros(1, size(X,2)); X];
Ylearn = [zeros(1, size(Y,2)); Y];


% ======================== variations
if doshift==1
    Xneural = Xneural(1:end-1,:);
    Ylearn = Ylearn(2:end,:);
%     Xneural = Xneural(2:end,:);
%     Ylearn = Ylearn(1:end-1,:);
elseif doshift==2
        Xneural = Xneural(2:end,:);
    Ylearn = Ylearn(1:end-1,:);
end

if bintobinchanges==1
    Xneural = diff(Xneural, 1, 1);
    Ylearn = diff(Ylearn, 1,1);
end



% ========== get unique learning trajectopries
[~, tmp] = intersect(exptID, exptID_U);
Ylearn_unique = Ylearn(:, tmp);

rhodat = diag(corr(Xneural, Ylearn))';
% lt_plot_histogram(rho_all);
% [~, p]= ttest(rho_all);
% % p= signrank(rho_all);
% lt_plot_pvalue(p, 'ttest',1);

% =============== SHUFFLE
rhoperm_all = nan(nshuff, length(rhodat));
for n=1:nshuff
    
    if (0)
        % ---- permute the order, keeping expt ID intact
        % PERMUTE ALL CASES
        Yperm = Ylearn(:, randperm(size(Ylearn,2)));
    else
        % PERMUTE AT LEVEL OF EXPERIMENT
        tmp = exptID_U(randperm(length(exptID_U)));
        Yperm = Ylearn_unique(:, tmp(exptID));
        %         Yperm = Ylearn_unique(:, exptID);
    end
    rhoperm = diag(corr(Xneural, Yperm));
    rhoperm_all(n, :) = rhoperm;
end


%% remove any nans
assert(all(find(isnan(rhodat)) == find(any(isnan(rhoperm_all)))), 'shuffle did not introduce errors (ie.. nans)');
indsnan = find(isnan(rhodat));
rhodat(indsnan) = [];
rhoperm_all(:, indsnan) = [];
