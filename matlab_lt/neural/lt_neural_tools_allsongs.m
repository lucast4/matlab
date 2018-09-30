function [Allbird_Fnames, Allbird_chanlist, Allbird_birdnum] = lt_neural_tools_allsongs(SummaryStruct)
%% lt 9/17/18 - extrast list of alls ongs and correspondings things (e.g. chans)


Allbird_Fnames = {};
Allbird_chanlist = {};
Allbird_birdnum = [];

for i=1:length(SummaryStruct.birds)
    disp(['bird ' num2str(i)]);
    
    birdname = SummaryStruct.birds(i).birdname;
    numneur = length(SummaryStruct.birds(i).neurons);
    
    All_Fnames = {};
    All_chans = [];
    for nn=1:numneur
        
        % ================== COLLECT LIST OF ALL FILENAMES FOR THIS NEURON
        dirname = SummaryStruct.birds(i).neurons(nn).dirname;
        indtmp = strfind(dirname, '/');
        dirname = dirname(1:indtmp(end)-1);
        datenames = SummaryStruct.birds(i).neurons(nn).Filedatestr_unsorted;
        
        fnamesall = {};
        for daten = datenames
            fnamethis = [dirname '/' birdname '_' daten{1} '.rhd'];
            fnamesall = [fnamesall; fnamethis];
        end
        
        chan = SummaryStruct.birds(i).neurons(nn).channel;
        chanvec = ones(size(fnamesall))*chan;
        
        % ==================== COLLECT FOR THIS NEUR
        All_Fnames = [All_Fnames; fnamesall];
        All_chans = [All_chans; chanvec];
    end
    
    
    % =========== GET MAPPING BETWEEN EACH SONG FILE AND LIST OF CHANS FOR IT
    FnamesUnique = unique(All_Fnames);
    ChanlistAll = {};
    for j=1:length(FnamesUnique)
        fnamethis = FnamesUnique{j};
        indthis = strcmp(All_Fnames, fnamethis);
        
        chanlistthis = All_chans(indthis);
        
        % ---- get unique chans (since sometimes a chan has SU and MU, or diff
        % clusts)
        assert(isempty(chanlistthis)==0, 'there must be chan...');
        chanlistthis = unique(chanlistthis);
        ChanlistAll = [ChanlistAll; chanlistthis];
    end
    
    % ######################## OUTPUT, ACROSS ALL BIRDS
    Allbird_Fnames = [Allbird_Fnames; FnamesUnique];
    Allbird_chanlist = [Allbird_chanlist; ChanlistAll];
    Allbird_birdnum = [Allbird_birdnum; ones(size(ChanlistAll))*i];

end

