function [Cohcell, Chanpairs, t, ff] = lt_neural_Coher_GetDat(fname_rhd, Chanlist)
% fname_rhd needs to be full dir (or in current dir)
% e.g.:
% fname_rhd = '/bluejay5/lucas/birds/pu69wh78/NEURAL/111317_RALMANOvernightLearn1/COHERENCE/pu69wh78_171113_173612.rhd'

%%  lt 10/1/18 - extract coherence for a given file, and list of channels

if ~exist([fname_rhd '.ffbins'], 'file')
    return
end

%%
% --- first sort chans so they are in order
Chanlist = sort(Chanlist);

% -------- GO THRU ALL
Cohcell = {};
Chanpairs = [];
for c=1:length(Chanlist)
    for cc=c+1:length(Chanlist)
        chan1 = Chanlist(c);
        chan2 = Chanlist(cc);
        
        % ---- coherence file
        fncoh = [fname_rhd '.ch' num2str(chan1) '-' num2str(chan2) '.coh'];
        datcoh = load(fncoh, '-mat');
        
        % ----- OUTPUT
        Cohcell = [Cohcell; datcoh.C];
        Chanpairs = [Chanpairs; [chan1 chan2]];
    end
end

% ---- load general things
fn_t = [fname_rhd '.tbins'];
fn_ff = [fname_rhd '.ffbins'];

t = load(fn_t, '-mat'); t = t.t;
ff = load(fn_ff, '-mat'); ff = ff.f;

% ============= OUTPUT STUFF;
