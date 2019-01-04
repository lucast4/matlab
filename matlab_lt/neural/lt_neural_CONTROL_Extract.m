function SummaryStruct = lt_neural_CONTROL_Extract(SummaryStruct)
%% lt 1/2/19 -first run summarrty struc.t then run this to keep only rthose designated by hand below

%% ################# DTABASE.

neurcontr = struct;
ind = 0;

% ==== here is database of channels to use [one ind for each directoruy];
ind = ind+1;
neurcontr(ind).birdname = 'pu69wh78';
neurcontr(ind).exptname = '110917_RALMANOvernightLearn1';
neurcontr(ind).chanstoget = [8 17];

ind = ind+1;
neurcontr(ind).birdname = 'pu69wh78';
neurcontr(ind).exptname = '110117_RALMANlearn1';
neurcontr(ind).chanstoget = [8 11 22];

ind = ind+1;
neurcontr(ind).birdname = 'pu69wh78';
neurcontr(ind).exptname = '110517_RALMANlearn2';
neurcontr(ind).chanstoget = [8 22];

ind = ind+1;
neurcontr(ind).birdname = 'pu69wh78';
neurcontr(ind).exptname = '103017_RAlearn1';
neurcontr(ind).chanstoget = [8 9 11 14 22];

% ----
ind = ind+1;
neurcontr(ind).birdname = 'wh44wh39';
neurcontr(ind).exptname = '021918_RALMANlearn1';
neurcontr(ind).chanstoget = [8 11 12 17 18 20 21 22];

ind = ind+1;
neurcontr(ind).birdname = 'wh44wh39';
neurcontr(ind).exptname = '031418_RALMANlearn2';
neurcontr(ind).chanstoget = [8 11 12 17 18 20 21 22];

ind = ind+1;
neurcontr(ind).birdname = 'wh44wh39';
neurcontr(ind).exptname = '032118_RALMANlearn3';
% neurcontr(ind).chanstoget = [8 11 12 17 18 20 21 22]; remove 22 since is
% short dur recortding (in  order to match actual dat).
neurcontr(ind).chanstoget = [8 11 12 17 18 20 21 22];

ind = ind+1;
neurcontr(ind).birdname = 'wh44wh39';
neurcontr(ind).exptname = '041718_RALMANlearn4';
neurcontr(ind).chanstoget = [8 11 12 17 18 20 21 22];

% ----
ind = ind+1;
neurcontr(ind).birdname = 'wh72pk12';
neurcontr(ind).exptname = '120518_RALMANLearn3';
neurcontr(ind).chanstoget = [10 18 23];

ind = ind+1;
neurcontr(ind).birdname = 'wh72pk12';
neurcontr(ind).exptname = '120718_RALMANLearn3';
neurcontr(ind).chanstoget = [10 13 18 23];

ind = ind+1;
neurcontr(ind).birdname = 'wh72pk12';
neurcontr(ind).exptname = '121018_RALMANLearn3';
neurcontr(ind).chanstoget = [10 13 18 23];

ind = ind+1;
neurcontr(ind).birdname = 'wh72pk12';
neurcontr(ind).exptname = '121718_RALMANLearn4';
neurcontr(ind).chanstoget = [10 13 23];

ind = ind+1;
neurcontr(ind).birdname = 'wh72pk12';
neurcontr(ind).exptname = '121918_RALMANLearn5';
neurcontr(ind).chanstoget = [10 13 18 23];


assert(length(lt_tools_grp2idx({{neurcontr.birdname}', {neurcontr.exptname}'})) == ind, 'problem entereing database, some not unique');

%% ================== DO FILTER (go thru all neurons. decide if keep or remove).
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    numneur = length(SummaryStruct.birds(i).neurons);
    neurtoremove = [];
    for ii=1:numneur
        dirthis = SummaryStruct.birds(i).neurons(ii).dirname;
        [~, enamethis] = fileparts(fileparts(dirthis));
%         disp(enamethis);
        chanthis = SummaryStruct.birds(i).neurons(ii).channel;
        birdthis = SummaryStruct.birds(i).birdname;
        
        % ==== throw out if is not wanted
        indtmp = strcmp({neurcontr.birdname}, birdthis) & strcmp({neurcontr.exptname}, enamethis);
            
        % ----- if not entered in database, then throw out
        if ~any(indtmp)
            neurtoremove = [neurtoremove, ii];
            continue
        else
            assert(sum(indtmp)==1, 'mistake in entering database...');
        end
        
        % ----- check if this channel wanted
        chanswanted = neurcontr(indtmp).chanstoget;
        if ismember(chanthis, chanswanted)
            % then keep
        else
            % then discard
            neurtoremove = [neurtoremove, ii];
        end
    end
    
    % ==== discard umnwanted neurons
    disp(['REMAINING NEURONS FOR ' birdthis]);
    disp([num2str(numneur - length(neurtoremove)) '/' num2str(numneur)]);
    SummaryStruct.birds(i).neurons(neurtoremove) = [];
        
end

%% ================== REEXTRACT POPULATION DATA
% ----- 1) REMOVE OLD DATA
SummaryStruct.birds =  rmfield(SummaryStruct.birds, 'exptnum_pop');

% ----- 2) GET POP[ DATA
SummaryStruct = lt_neural_PRE_GetSimultNeur(SummaryStruct);



