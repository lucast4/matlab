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

%% ================ MAKE BATCH FILES WITH LABELED SONGS.
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
overwrite = 0; % if 0, then will do any cases not already done.
plotEachSyl = 0;

for i=1:nbirds
    ndir = length(dirstruct.bird(i).DIR_directories);
    for ii=1:ndir
        
        % =============================================
        ListOfDirs_UNDIR = {dirstruct.bird(i).UNDIR_directories{ii}};
        ListOfDirs_DIR = {dirstruct.bird(i).DIR_directories{ii}};
        
        ListOfDirs_ALL = [ListOfDirs_UNDIR ListOfDirs_DIR];
        
        ListOfBatch = {...
            'batchall.LABELED', ...
            'batchall.LABELED'};
        
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
DATSTRUCT = lt_batchsong_extractFF(ListOfDirs_UNDIR, ListOfDirs_DIR, ...
    ListOfBatch, MotifsToExtract);



%% ===================== FOR EACH BIRD, EXTRACT ALL UNIQUE MOTIFS AND SAVE
dirstruct = lt_DirSong_MotifID;


%% ================ PLOT SUMMARY FOR EACH BIRD

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

nbirds = length(dirstruct.bird);
for i=1:nbirds
    nmotifs = length(dirstruct.bird(i).DAT.motifID);
    for mm=1:nmotifs
        
        disp([i mm]);
        ff = [dirstruct.bird(i).DAT.motifID(mm).rendnum.ff];
        t = [dirstruct.bird(i).DAT.motifID(mm).rendnum.datenum_song_SecRes];
        isDir = [dirstruct.bird(i).DAT.motifID(mm).rendnum.isDIR];
        
        if all(isnan(ff))
            disp('ff is nan...');
            continue
        end
        
        % === sort by t (for fun)
        [~, indsort] = sort(t);
        t = t(indsort);
        ff = ff(indsort);
        isDir = isDir(indsort);
        
        if rand<1.1
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([bname '-mot' num2str(mm)]);
            plot(t(isDir==0), ff(isDir==0), 'ok');
            plot(t(isDir==1), ff(isDir==1), 'ob');
            datetick('x', 'mm/dd');
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([bname '-mot' num2str(mm)]);
            
            x = 1:length(t);
            plot(x(isDir==0), ff(isDir==0), 'ok');
            plot(x(isDir==1), ff(isDir==1), 'ob');
        end
        
        
        % ========== for each day, get mean DIR and UNDIR (difference)
        t_day = floor(t);
        
        [ff_day_UNDIR, ff_day_UNDIR_std] = grpstats(ff(isDir==0), t_day(isDir==0), {'mean', 'std'});
        [ff_day_DIR, ff_day_DIR_std] = grpstats(ff(isDir==1), t_day(isDir==1), {'mean', 'std'});
        %     ff_day_UNDIR = grpstats(ff(isDir==0), t_day(isDir==0), {'median'});
        %     ff_day_DIR = grpstats(ff(isDir==1), t_day(isDir==1), {'median'});
        assert(length(ff_day_UNDIR)==length(ff_day_DIR));
        
        ff_UndirMinusDir = ff_day_UNDIR - ff_day_DIR;
        
        %     ff_UndirMinusDir_z = [];
        ff_UndirMinusDir_z = -(ff_day_DIR - ff_day_UNDIR)./ff_day_UNDIR_std;
        
        
        % ============= COULD BE GOOD:
%         % ==== CV UNDIR AND DIR
%         cv_UNDIR = ff_day_UNDIR_std./ff_day_UNDIR;
%         cv_DIR = ff_day_DIR_std./ff_day_DIR;
%         
%         All_ffUndir_cv = [All_ffUndir_cv; cv_UNDIR];
%         All_ffDir_cv = [All_ffDir_cv; cv_DIR];
%         
%         if strcmp(cvdiffmethod, 'diff')
%             cvdiff_UndirOverDir = cv_UNDIR-cv_DIR;
%         elseif strcmp(cvdiffmethod, 'percent')
%             cvdiff_UndirOverDir = (cv_UNDIR-cv_DIR)./cv_UNDIR;
%         end
%         All_ffUndirOverDir_cv = [All_ffUndirOverDir_cv; cvdiff_UndirOverDir];
%         
%         % =============== SAVE
%         if useff_zscore==2
%             All_ffUndirMinusDir = [All_ffUndirMinusDir; (ff_day_UNDIR - ff_day_DIR)./ff_day_UNDIR];
%         elseif useff_zscore==1
%             All_ffUndirMinusDir = [All_ffUndirMinusDir; ff_UndirMinusDir_z];
%         elseif useff_zscore==0
%             All_ffUndirMinusDir = [All_ffUndirMinusDir; ff_UndirMinusDir];
%         end
%         All_neuralHiMinusLo = [All_neuralHiMinusLo; yneur];
%         All_bnum = [All_bnum; bnum];
%         All_motifID = [All_motifID; mm];
%         
%         All_neuralXcovBase = [All_neuralXcovBase; yxcov_base];
    end
end