function SummaryStructOut = lt_neural_RASamMel_SummaryStruct

%% ======== make Summary struct
cd('/bluejay5/lucas/Sober_Mel_RA_physiology/RA_data_summaries_putative_proj_neurons');

datfiles = dir('combined_data*.mat');


%% ======================= 1) get list of birds
birdnamesAll = {};
for i=1:length(datfiles)
    uscores = strfind(datfiles(i).name, '_');
    birdname = datfiles(i).name(uscores(2)+1:uscores(3)-1);
    
    birdnamesAll = [birdnamesAll birdname];
end

BirdList = unique(birdnamesAll);


%% ====================== 2) slide into Summarystruct
SummaryStructOut = struct;
AllUseUnit = [];

for i=1:length(BirdList)
    birdname = BirdList{i};
    
    % ------------- go thru all dat files. skip if not this bird
    neuroncount = 0;
    for ii=1:length(datfiles)
        uscores = strfind(datfiles(ii).name, '_');
        birdnametmp = datfiles(ii).name(uscores(2)+1:uscores(3)-1);
        
        
        if strcmp(birdname, birdnametmp)==0
            continue
        end
        
        % ========== collect date
        tmp1 = strfind(datfiles(ii).name, '_MU_');
        tmp2 = strfind(datfiles(ii).name, '_SU_');
        
        if ~isempty([tmp1 tmp2])
            % then name version 1
            dateofneur = datfiles(ii).name(uscores(4)+1:uscores(7)-1);
        else
            dateofneur = datfiles(ii).name(uscores(3)+1:uscores(6)-1);
        end
        disp(dateofneur);
        disp(datfiles(ii).name);
        
        % ============= COLLECT DAT
        neuroncount = neuroncount+1;
        
        fname = [pwd '/' datfiles(ii).name];
        
        tmp = load(fname);
        assert(strcmp(tmp.birdname, birdname)==1, 'asdf');
        
        % =================== COLLECTING THINGS FOR MY CURIOSITY
        AllUseUnit = [AllUseUnit tmp.use_unit];
        
        % ============= any ff that is 0 - convert to nan
        tmp.peak_pinterp_labelvec(tmp.peak_pinterp_labelvec==0)=nan;
        
       
        
        % ============= OUTPUT
        SummaryStructOut.birds(i).birdname = birdname;
        SummaryStructOut.birds(i).neurons(neuroncount).datfilename = fname;
        SummaryStructOut.birds(i).neurons(neuroncount).date = dateofneur;
        SummaryStructOut.birds(i).neurons(neuroncount).dat = tmp;
        
    end
    
end

% %% =============
% lt_figure; hold on;
