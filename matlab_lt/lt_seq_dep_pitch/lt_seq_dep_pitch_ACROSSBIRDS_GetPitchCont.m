function lt_seq_dep_pitch_ACROSSBIRDS_GetPitchCont(SeqDepPitch_AcrossBirds, ...
    birdtoget, expttoget)
%% lt 7/10/18 - extracts all pitch contours, aligns to trials

%%

fdirmain = '/bluejay5/lucas/analyses/seqdeppitch/GetPitchCont';

%% ============ go to folder with raw dat

i = birdtoget;
ii = expttoget;

fname = [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/AllDays_RawDatStruct.mat'];
fname_params = [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/Params.mat'];

load(fname);
load(fname_params);


%% ============= from rawdat keep only needed things (fre up memry)


%% ============== first make directory to save for this bird and expt
birdname = SeqDepPitch_AcrossBirds.birds{i}.birdname;
exptname = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;

disp([birdname '-' exptname]);

fdirthis = [fdirmain '/' birdname '_' exptname];

if exist(fdirthis, 'dir')==0
    mkdir(fdirthis);
end

%% ============= EXTRACT ALL PC

SylList = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
Numsyls = length(SylList);

for ss = 1:Numsyls
    
    sylthis = SylList{ss};
    Datthis = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii};
    Numdays = length(Datthis.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis).FFvals);
    
    assert(Numdays == length(AllDays_RawDatStruct));
    
    for dd = 1:Numdays
        
        % ------- is today empty?
        if isempty(AllDays_RawDatStruct{dd})
            continue
        end
        
        PCmat = cell2mat(AllDays_RawDatStruct{dd}.data.(sylthis)(:,2));
        FFmat = cell2mat((AllDays_RawDatStruct{dd}.data.(sylthis)(:,1)));
        Tmat = cell2mat((AllDays_RawDatStruct{dd}.data.(sylthis)(:,6)));
        
            PCmat = [cell2mat(AllDays_RawDatStruct{dd}.data_MUSC.(sylthis)(:,2)); ...
                cell2mat(AllDays_RawDatStruct{dd}.data.(sylthis)(:,2))];
            
            FFmat = [cell2mat(AllDays_RawDatStruct{dd}.data_MUSC.(sylthis)(:,1)); ...
                cell2mat(AllDays_RawDatStruct{dd}.data.(sylthis)(:,1))];
            Tmat = [cell2mat(AllDays_RawDatStruct{dd}.data_MUSC.(sylthis)(:,6)); ...
                cell2mat(AllDays_RawDatStruct{dd}.data.(sylthis)(:,6))];

            % -------- GET INDS, compare to previously extracted data
        FFold = Datthis.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis).FFvals_WithinTimeWindow{dd};
        Told = Datthis.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis).Tvals_WithinTimeWindow{dd};
        %         Datthis.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(sylthis).FFvals{dd}
%         [~, ~, indtmp] = intersect(FFold+Told, FFmat+Tmat, 'stable');
        indtmp = find(ismember(Tmat, Told));
        assert(length(unique(indtmp))==length(FFold)); assert(all(diff(indtmp)>0));
        assert(all(Told' == Tmat(indtmp)));
        assert(all(FFold' == FFmat(indtmp)))
        
        % ------------- to save
        PCmat_tosave = single(PCmat(indtmp,:));
        fname = [fdirthis '/syl' sylthis '_PBS_day' num2str(dd) '.mat'];
        save(fname, 'PCmat_tosave');
        
        
        
        % ================== COLLECT MUSCIMOL DAT?
        FFold_MUSC = Datthis.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(sylthis).FFvals_WithinTimeWindow{dd};
        Told = Datthis.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(sylthis).Tvals_WithinTimeWindow{dd};
        if ~isempty(FFold_MUSC)
            
            % === NOTE: include all trials for the day, since some MUSC
            % data comes at transition between MUSc and PBS.
            PCmat_MUSC = [cell2mat(AllDays_RawDatStruct{dd}.data_MUSC.(sylthis)(:,2)); ...
                cell2mat(AllDays_RawDatStruct{dd}.data.(sylthis)(:,2))];
            
            FFmat_MUSC = [cell2mat(AllDays_RawDatStruct{dd}.data_MUSC.(sylthis)(:,1)); ...
                cell2mat(AllDays_RawDatStruct{dd}.data.(sylthis)(:,1))];
            Tmat_MUSC = [cell2mat(AllDays_RawDatStruct{dd}.data_MUSC.(sylthis)(:,6)); ...
                cell2mat(AllDays_RawDatStruct{dd}.data.(sylthis)(:,6))];
            
            % --------------------- align renditions to old data
%             [~, ~, indtmp] = intersect(FFold_MUSC+Told, FFmat_MUSC+Tmat_MUSC, 'stable');
            indtmp = find(ismember(Tmat_MUSC, Told));
            assert(length(unique(indtmp))==length(FFold_MUSC)); assert(all(diff(indtmp)>0));
            assert(all(Told' == Tmat_MUSC(indtmp)));
            assert(all(FFold_MUSC' == FFmat_MUSC(indtmp)))
            
            % --- save
            PCmat_tosave_MUSC = single(PCmat_MUSC(indtmp,:));
            fname = [fdirthis '/syl' sylthis '_MUSC_day' num2str(dd) '.mat'];
            save(fname, 'PCmat_tosave_MUSC');
        end
    end
end


if (0)
    %% ad hoc, check a specific day and syl. modify this code for general code.
    
    
    PCall = cell2mat(AllDays_RawDatStruct{4}.data.ljbB(:,2));
    FFall = cell2mat(AllDays_RawDatStruct{4}.data.ljbB(:,1));
    
    
    lt_figure; hold on;
    lt_plot_histogram(FFall);
    
    lt_figure; hold on;
    
    % ---- sort in order
    [~, indstmp] = sort(FFall);
    FFall = FFall(indstmp);
    PCall = PCall(indstmp,:);
    
    N=24;
    xwind = 20:180;
    % --------------- lowest FF
    lt_subplot(2,2,1); hold on;
    title('low FF');
    pcthis = PCall(1:N,:);
    ffthis = FFall(1:N);
    % - subtract mean across trials
    tmp = mean(pcthis,1);
    pcthis = pcthis - repmat(tmp, N, 1);
    % - subtract mean within trial
    tmp = mean(pcthis(:,xwind),2);
    pcthis = pcthis - repmat(tmp,1,size(pcthis,2));
    
    plot(pcthis', '-k');
    
    % ---- get std within window
    lt_subplot(2,2,2); hold on;
    ylabel('std');
    y = std(pcthis(:, xwind), [], 2);
    % norm to ff
    y = y./ffthis;
    lt_plot_histogram(y);
    y1 = y;
    
    % --------------- lowest FF
    lt_subplot(2,2,3); hold on;
    title('high FF');
    pcthis = PCall(end-N+1:end, :);
    ffthis = FFall(end-N+1:end);
    % - subtract mean across trials
    tmp = mean(pcthis,1);
    pcthis = pcthis - repmat(tmp, N, 1);
    % - subtract mean within trial
    tmp = mean(pcthis(:,xwind),2);
    pcthis = pcthis - repmat(tmp,1,size(pcthis,2));
    
    plot(pcthis', '-k');
    
    % ---- get std within window
    lt_subplot(2,2,4); hold on;
    ylabel('std');
    y = std(pcthis(:, xwind), [], 2);
    % norm to ff
    y = y./ffthis;
    lt_plot_histogram(y);
    
    %
    % % ######################### USING cv
    % lt_figure; hold on;
    % % --------------- lowest FF
    % lt_subplot(2,2,1); hold on;
    % title('low FF');
    % pcthis = PCall(1:N,:);
    %
    % % -------- take time of interest
    % pcthis = pcthis(:,xwind);
    % % ------- subtract mean for each trial
    %
    %
    % % - subtract mean across trials
    % tmp = mean(pcthis,1);
    % pcthis = pcthis - repmat(tmp, N, 1);
    % % - subtract mean within trial
    % tmp = mean(pcthis(:,xwind),2);
    % pcthis = pcthis - repmat(tmp,1,size(pcthis,2));
    %
    % plot(pcthis', '-k');
    
    %% OLD, IN PROGRESS.
    
    sylthis = 'jjB';
    daythis = 1;
    
    datmat = tmp.AllDays_RawDatStruct{daythis}.data.(sylthis);
    
    figure; hold on;
    plot(datmat{1,2}, '-');
    
    % ==== get window for contour...
    
    
    % ==== get contours for all trials, aligned to the trials used in the
    % analysis
    
    
    Datthis = SeqDepPitch_AcrossBirds.birds{10}.experiment{1}.Data_PlotLearning;
    Datthis.AllDays_PlotLearning.DataMatrix.(sylthis).FFvals
end
