function F = lt_tools_lme_plotPredict(lme, X)
%% plot predictions conditioned on training data (X empyt) or new data (X)

if ~exist('X', 'var')
    X = [];
end

%% 
if isempty(X)
    F = lme.fitted;
    X = lme.Variables;
else
    F = lme.predict(X);
end

%% === 

lme.CoefficientNames

