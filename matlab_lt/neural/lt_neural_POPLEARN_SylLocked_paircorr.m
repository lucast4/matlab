function rhoall = lt_neural_POPLEARN_SylLocked_paircorr(frmat)
%%  2/26/19 - LT, given frmat outputs all pairwise correlations
% frmat is cell array, with each cell diff unit. within each cell it is
% timebin x trial

% - output:
% rhoall is trial x nejron paur;


% disp('NOTE: any cases with nan for correlation I replace by median acros tirlas (for that neuron pair x syllable)');

%%

% ========== go thru all trials
ntrials = size(frmat{1},2);
nunits = length(frmat);
rhoall = nan(ntrials, nchoosek(nunits, 2)); % trials x pair
for j=1:ntrials
    % -- go thru all spiking pairs
    cc = 1;
    for nn=1:nunits
        for nnn=nn+1:nunits
            fr1 = frmat{nn}(:,j);
            fr2 = frmat{nnn}(:,j);
            
            % -- get cross correlation
            rhoall(j, cc) = corr(fr1, fr2);
            cc = cc+1;
        end
    end
end
% ---- check for all nan
if sum(isnan(rhoall(:)))/length(isnan(rhoall(:)))>0.01
    disp('more than 1% cases got nan for correlation ...');
    keyboard
end
% --- convert nan to median
ncols = size(rhoall,2);
for j=1:ncols
    rhoall(isnan(rhoall(:,j)), j) = nanmedian(rhoall(:,j));
end

