%% Directed song database
% --- extraction of FF, and database of birds.



%% ------ LIST DIRECTORIES OF ALL DIRECTED SONG
clear all; close all;
dirstruct = lt_DirSong_Extract;




%% ================ run in a given folder, do to check what is good params
if (0)
    close all;
    ListOfDirs_UNDIR = {...
        '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/UNDIR'};
    
    
    ListOfDirs_DIR = {...
        };
    
    ListOfDirs_ALL = [ListOfDirs_UNDIR ListOfDirs_DIR];
    
    ListOfBatch = {...
        'batchall'};
    
    FFparams.cell_of_freqwinds={'a', [700 2600], 'g', [2000 3400], ...
        'h', [2900 3800], 'b', [2900 3800]}; % 'j', [950 1450], 'l', [1200 1600], 't', [3590 4960]
    FFparams.cell_of_FFtimebins={'a', [0.068 0.078], 'g', [0.04 0.06], ...
        'h', [0.033 0.042], 'b', [0.033 0.04]}; % 'j', [0.04 0.045], 'l', [0.035 0.039], 't', [0.026 0.033], ...
    
    plotAllPC = 1;
    plotEachSyl = 0;
    overwrite = 1;
    
    
    % ========= generate batch
    cd(ListOfDirs_UNDIR{1});
    eval('!ls *.rhd > batchall');
    eval('!ls *.cbin >> batchall');
    
    % ==================== CALCULATE AND SAVE FF
    lt_batchsong_calcFF(ListOfDirs_ALL, ListOfBatch, FFparams, plotAllPC, plotEachSyl, ...
        overwrite);
end

%% ================ EXTRACTION OF FF FROM SONG FILES
nbirds = length(dirstruct.bird);

% === will skip if already done

% =============== 1) in each folder make a batch of all the songs.
for i=1:nbirds
    ndir = length(dirstruct.bird(i).DIR_directories);
    for ii=1:ndir
        % ================ DIRECTED
        cd(dirstruct.bird(i).DIR_directories{ii});
%         disp(dirstruct.bird(i).DIR_directories{ii});
        eval('!> batchall');
        eval('!ls *.rhd > batchall');
        eval('!ls *.cbin >> batchall');
% ----- GET SUBNSET THAT IS LABELED
lt_sort_batch_by_labeled('batchall');
        
        
        cd(dirstruct.bird(i).UNDIR_directories{ii});
        eval('!> batchall');
        eval('!ls *.rhd > batchall');
        eval('!ls *.cbin >> batchall');
% ----- GET SUBNSET THAT IS LABELED
lt_sort_batch_by_labeled('batchall');
        

    end
end


%% =============== EXTRACT FF FOR EACH CASE
nbirds = length(dirstruct.bird);

% doOverwrite = 1;
plotAllPC = 1;
overwrite = 1; % if 0, then will do any cases not already done.
plotEachSyl = 0;

for i=1:nbirds
    ndir = length(dirstruct.bird(i).DIR_directories);
    for ii=1:ndir
        
        % =============================================
        ListOfDirs_UNDIR = {dirstruct.bird(i).UNDIR_directories{ii}};
        ListOfDirs_DIR = {dirstruct.bird(i).DIR_directories{ii}};
        
        ListOfDirs_ALL = [ListOfDirs_UNDIR ListOfDirs_DIR];
        
        ListOfBatch = {...
            'batchall.UNLABELED', ...
            'batchall.UNLABELED'};
        
        FFparams = dirstruct.bird(i).FFparams;
        
        
        % ========= generate batch
        %         cd(ListOfDirs_UNDIR{1});
        %         eval('!ls *.rhd > batchall');
        %         eval('!ls *.cbin >> batchall');
        %
        % ==================== CALCULATE AND SAVE FF
        lt_batchsong_calcFF(ListOfDirs_ALL, ListOfBatch, FFparams, plotAllPC, plotEachSyl, ...
            overwrite);
        lt_save_all_figs;
        close all;
    end
end



%% ###########################################################
%% ###########################################################
%% ==================== EXTRACT FF FOR A GIVEN SYLLABLE
MotifsToExtract = {'ab(h)', 'jb(h)',  'jbh(h)', 'a(b)h', 'j(b)', 'h(g)'};
DATSTRUCT = lt_batchsong_extractFF(ListOfDirs_UNDIR, ListOfDirs_DIR, ListOfBatch, MotifsToExtract);



%% ===================== FOR EACH BIRD, EXTRACT ALL UNIQUE MOTIFS AND SAVE
dirstruct = lt_DirSong_MotifID;