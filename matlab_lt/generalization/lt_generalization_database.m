%% lt 4/4/18 - stores information about each experiment.
function [DATSTRUCT, TargSyl, MotifList] = lt_generalization_database(birdname, exptname)

%%
if strcmp(birdname, 'pu69wh78') & strcmp(exptname, 'RALMANlearn1')
    
    % --- PARAMS
    ListOfDirs_ALL = {...
        '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1'};
    ListOfBatch = {...
        'BatchAll.LABELED'};
    MotifsToExtract = {'a(b)', 'j(b)', 'ab(h)', 'jb(h)',  'jbh(h)', '(g)'};
    TargSyl = 'a(b)';
    % --- EXTRACT
    DATSTRUCT = lt_batchsong_extractFF(ListOfDirs_ALL, {}, ListOfBatch, MotifsToExtract);
    
end


%% 
MotifList = MotifsToExtract;

%% seq dep pitch


