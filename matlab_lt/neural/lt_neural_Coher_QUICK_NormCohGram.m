function cohmatNorm = lt_neural_Coher_QUICK_NormCohGram(cohmatAll)
% ==== TO LOOK AT FLUCTIOATOINS, NORMS COHGRAM 
% ---- SBTRACTS ACROSS TIME (WITHIN F) MEAN, SEPARATELY FOR EACH IND (DIM
% 3)

% cohmat = t x ff x trial(or syl, or whatever)

%% 

cohmatNorm = nan(size(cohmatAll));
for i=1:size(cohmatAll,3)
    
    cohmat = cohmatAll(:,:, i);
    cohmeanthis = mean(cohmat,1);
    cohmat = cohmat - repmat(cohmeanthis, size(cohmatAll,1), 1);
    
    cohmatNorm(:,:,i) = cohmat;    
end
