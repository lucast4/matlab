function lt_batchsong_plotallchans(batchf, lasernote, annotate_note_group)
%% lt 9/27/18 - modified from previuos code used specificalyl for opto experiemnts 
% plots song, other chans (e.g laser) and times of hits.

% batchf, songfile or batch
% lasernote = 0, 1, etc
% annotate_note_group = 1; % if 1, then shows which NG each song was in.

songf = batchf; % conversion from old code.

%%
if any(strfind(songf, 'cbin'))
    lt_figure; hold on;
    % then assume is a single song
    % 1) song
    [dat, Fs, DOFILT, ext]=ReadDataFile(songf,'0');
    x = [1:length(dat)]./Fs;
    plot(x, dat, '-k')
    
    % 2) trigger
    [dat, Fs, DOFILT, ext]=ReadDataFile(songf,'1');
    plot(x, dat, '-r')
    
    % -- ad trig time with rec file ...
    rd = readrecf_LT_evtafv4(songf);
    
    
else
    
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    
    hsplots = [];
    % ---------- go thur all songs in batch
    batchf = songf;
    fid = fopen(batchf);
    songf = fgetl(fid);
    
    while ischar(songf)
        
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(songf);
        hsplots = [hsplots hsplot];
        
        % 1) song
        [dat, Fs, DOFILT, ext]=ReadDataFile(songf,'0');
        x = [1:length(dat)]./Fs;
        plot(x, dat, '-k')
        
        % 2) trigger
        [dat, Fs, DOFILT, ext]=ReadDataFile(songf,'1');
        plot(x, dat, '-r')
        
        % -- ad trig time with rec file ...
        rd = readrecf_LT_evtafv4(songf);
        if ~isempty(rd.ttimes)
            for i=1:length(rd.ttimes)
                if rd.trignote(i)==lasernote
                    if rd.catch(i) ==1
                        plot(rd.ttimes(i)/1000, 0, '^g');
                    else
                        plot(rd.ttimes(i)/1000, 0, '^r');
                        lt_plot_text(rd.ttimes(i)/1000 +0.05, 0, 'hit', 'r');
                    end
                    
                end
                
            end
        end
        
        % - next
        songf = fgetl(fid);
    end
    linkaxes(hsplots, 'y');
end