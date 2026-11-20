function run_xcorr_shape(deployment,outDisk,fileMask)
% RUN_XCORR_SHAPE  Run the shape-based cross-correlation detector on a whole deployment.
%
%   run_xcorr_shape('GOM_MC_17','E:\AirgunShape')
%   run_xcorr_shape('\\frosty.ucsd.edu\GofMX_Decimated_12\GOM_GC_18','E:\AirgunShape')
%
% deployment: folder name under \\frosty.ucsd.edu\GofMX_Decimated_8, e.g. 'GOM_MC_17',
%             or the full path of a deployment folder on any share.
%             Every df100 xwav in every disk folder of that deployment is processed.
% outDisk:    where results go. One folder per deployment is created,
%             <outDisk>\<deployment>_AirgunDetections, with one .mat per xwav
%             (same variables as xcorr_airgun_df100, plus allTier and allCand).
% fileMask:   optional; process only files matching this pattern (for a quick test).
%
% Files that already have a result are skipped, so a run that is interrupted
% can simply be started again. To run several deployments at once, open one
% MATLAB session per deployment.
%
% Needs Triton's utFindFiles on the path (Triton\Remoras\Logger).

here = fileparts(mfilename('fullpath')); repo = fileparts(here);
root = '\\frosty.ucsd.edu\GofMX_Decimated_8';
if any(deployment == '/' | deployment == '\')  % a full path was given
    [root,name,ext] = fileparts(deployment); deployment = [name ext];
end
addpath(fullfile(repo,'xcorr_original'));   % readxwavhd and the template live with the original detector
if isempty(which('utFindFiles'))
    addpath(fullfile(fileparts(repo),'Triton','Remoras','Logger'));
end
assert(~isempty(which('utFindFiles')),'utFindFiles not found: add Triton\\Remoras\\Logger to the path.');
assert(exist(fullfile(root,deployment),'dir') == 7,'Deployment folder not found: %s',fullfile(root,deployment));
if ~exist(outDisk,'dir'), mkdir(outDisk); end
if nargin < 3, fileMask = '*df100.x.wav'; end

% a runnable copy of the prototype, changing only where it reads and writes
src = fileread(fullfile(here,'xcorr_airgun_shape.m'));
sub = {
    '^function\s+\w+\s*$',        'function ag_prototype_run'
    '^deplList = [^\n]*$',         sprintf('deplList = dir(''%s'');',fullfile(root,[deployment '*']))
    '^outDisk\s*= [^\n]*$',        sprintf('outDisk = ''%s'';',outDisk)
    '^templateFile = [^\n]*$',     sprintf('templateFile = ''%s'';',fullfile(repo,'xcorr_original','air_template_df100.mat'))
    '^SearchFileMask = [^\n]*$',   sprintf('SearchFileMask = {''%s''};',fileMask)
    '^overWritePrevious = [^\n]*$','overWritePrevious = 0;'};
for i = 1:size(sub,1)
    assert(numel(regexp(src,sub{i,1},'lineanchors','dotexceptnewline')) == 1,'could not find the line to adapt: %s',sub{i,1});
    src = regexprep(src,sub{i,1},strrep(sub{i,2},'\','\\'),'lineanchors','dotexceptnewline');
end
work = tempname; mkdir(work); addpath(work);
fid = fopen(fullfile(work,'ag_prototype_run.m'),'w'); fwrite(fid,src); fclose(fid);
cleanup = onCleanup(@() rmpath(work));

fprintf('Shape-based detector: %s -> %s\n',deployment,fullfile(outDisk,[deployment '_AirgunDetections']));
t0 = tic; ag_prototype_run; fprintf('Finished %s in %.1f hours.\n',deployment,toc(t0)/3600);
end
