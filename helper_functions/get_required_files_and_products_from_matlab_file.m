%get required files and products from matlab file


%% for PFC paper (2021)

[fList,pList] = matlab.codetools.requiredFilesAndProducts('RunAnalysis.m');
newplace = 'C:\Users\Alice Berners-Lee\Dropbox\FosterLab\scripts\For_RunAnalysis\';
if ~isfolder(newplace)
    mkdir(newplace)
end
for idir = 1:size(fList,2)
    oldplace = fList{idir};
    slashes = strfind(oldplace,'\');
    thename = oldplace(slashes(end)+1:end);
    
    copyfile(oldplace,[newplace thename])
end

%% for Target paper (2022)

%%% need to add path
addpath(genpath('C:\Users\AliceBL\Dropbox\FosterLab\scripts'))

[fList,pList] = matlab.codetools.requiredFilesAndProducts(...
    'C:\Users\AliceBL\Dropbox\FosterLab\scripts\MurthyLab\LizData\Paper\TargetPaper_ReadMe_AllAnalysis_AllFigures.m');
newplace = 'C:\Users\AliceBL\Dropbox\MurthyLab\TargetPaper\Scripts_For_TargetPaper_edits\';


%%
if ~isfolder(newplace)
    mkdir(newplace)
end
for idir = 1:size(fList,2)
    oldplace = fList{idir};
    slashes = strfind(oldplace,'\');
    thename = oldplace(slashes(end)+1:end);
    
    copyfile(oldplace,[newplace thename])
end
