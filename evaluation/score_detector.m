function T = score_detector(outName,overrides,detectorFile,setName)
% SCORE_DETECTOR  Run a detector on the windows of a labelled set and score it.
%
% Runs the detector exactly as written in ../xcorr_original/xcorr_airgun_df100.m (same
% settings and logic), limited to the four segments of each test window,
% and compares its detections with the reviewed shot times.
%
%   found   shots that fall inside a detection (within tol seconds)
%   missed  shots with no detection
%   merged  shots that share a detection with an earlier shot (counted in missed)
%   false   detections that contain no shot
%
% Truth: windows whose label says "no airguns" have zero shots; for the
% others the reviewed shot times in truth.csv are used (one event per shot
% cycle, even when a cycle has several overlapping arrivals).
%
% outName:   name for the results file, e.g. 'baseline' -> scores_baseline.csv
% overrides: optional N-by-2 cell of {setting, value-as-text} to try other
%            settings without editing the detector, e.g.
%            {'minCorrThresh','1e15'; 'parm.c2_offset','1'}
% detectorFile: optional path to another detector to score (default
%            ../xcorr_original/xcorr_airgun_df100.m)
% setName:   folder of the labelled set, e.g. 'set2' (default 'set1')
%
% When the set has a review.csv (written by review_marks), it decides what
% each window is: "no airguns" windows have zero shots, "dropped" windows
% are left out, and the rest use the marks in truth.csv. "unknown" windows
% (the reviewer was not sure) are scored as no airguns and keep that label
% in the status column, so they can be reported as a caveat.

if nargin < 1, outName = 'baseline'; end
if nargin < 2, overrides = cell(0,2); end
if nargin < 4 || isempty(setName), setName = 'set1'; end
if nargin < 3 || isempty(detectorFile), detectorFile = fullfile(fileparts(fileparts(mfilename('fullpath'))),'xcorr_original','xcorr_airgun_df100.m'); end
fs = 2000; N = 150000; blk = 4; tol = 2;      % tol: seconds of slack when matching a shot to a detection
here = fileparts(mfilename('fullpath')); repo = fileparts(here); setDir = fullfile(here,setName);
root = '//frosty.ucsd.edu/GofMX_Decimated_8';
addpath(fullfile(repo,'xcorr_original'));
if isempty(which('utFindFiles'))               % the detector needs Triton's file finder
    addpath(fullfile(fileparts(repo),'Triton','Remoras','Logger'));
end

% ---- make a runnable copy of the detector, changing only where it reads and writes ----
src = fileread(detectorFile);
sub = {
    '^function\s+\w+\s*$',                'function ag_under_test'
    '^deplList = [^\n]*$',                 'deplList = dir(getenv(''AG_DEPL''));'
    '^outDisk\s*= [^\n]*$',                'outDisk = getenv(''AG_OUT'');'
    '^templateFile = [^\n]*$',             'templateFile = getenv(''AG_TEMPLATE'');'
    '^SearchFileMask = [^\n]*$',           'SearchFileMask = {getenv(''AG_MASK'')};'
    '^(\s*)for idx = 1:inc\s*$',           '$1for idx = str2double(getenv(''AG_SEG0'')):str2double(getenv(''AG_SEG1''))'
    '^\s*display\(\[''(start of segment|calculating)[^\n]*$', ''};
for i = 1:size(sub,1)
    n = numel(regexp(src,sub{i,1},'lineanchors','dotexceptnewline'));
    assert(n >= 1,'could not find the line to adapt: %s',sub{i,1});
    src = regexprep(src,sub{i,1},sub{i,2},'lineanchors','dotexceptnewline');
end
for i = 1:size(overrides,1)
    pat = ['^(\s*)' regexptranslate('escape',overrides{i,1}) '\s*=\s*[^;\n]*;'];
    assert(numel(regexp(src,pat,'lineanchors','dotexceptnewline')) == 1,'setting not found exactly once: %s',overrides{i,1});
    src = regexprep(src,pat,['$1' overrides{i,1} ' = ' overrides{i,2} ';'],'lineanchors','dotexceptnewline');
end
work = tempname; mkdir(work); addpath(work);
fid = fopen(fullfile(work,'ag_under_test.m'),'w'); fwrite(fid,src); fclose(fid);
cleanup = onCleanup(@() rmpath(work));
setenv('AG_TEMPLATE',fullfile(repo,'xcorr_original','air_template_df100.mat'));

W = readtable(fullfile(setDir,'windows.csv'),'Delimiter',',','TextType','string');
P = readtable(fullfile(setDir,'truth.csv'),'Delimiter',',','TextType','string');
if ~isstring(P.window), P.window = string(P.window); end      % an empty truth file is read without text columns
if exist(fullfile(setDir,'review.csv'),'file')
    Rv = readtable(fullfile(setDir,'review.csv'),'Delimiter',',','TextType','string');
    [~,ix] = ismember(W.window,Rv.window); W.status = Rv.status(ix);
else
    W.status = repmat("airguns",height(W),1); W.status(contains(W.condition_guess,'no airguns')) = "no airguns";
end
W = W(~ismember(W.status,["dropped","airguns unclear"]),:);   % "airguns unclear": airguns present but shots could not be marked one by one
% humps.csv (optional) lists windows where the reviewer could not tell whether the pulses in a
% hump are separate shots or arrivals of one shot: marks closer than join_s count as one event
Hp = table(strings(0,1),zeros(0,1),'VariableNames',{'window','join_s'});
if exist(fullfile(setDir,'humps.csv'),'file')
    Hp = readtable(fullfile(setDir,'humps.csv'),'Delimiter',',','TextType','string');
end   % "unknown" windows stay in, scored as no airguns, and are flagged in the status column
nW = height(W); res = zeros(nW,7); mins = cell(0,4); dets = cell(0,3); feats = cell(0,11); cands = cell(0,11);
for i = 1:nW
    d = dir(fullfile(root,"GOM_" + W.site(i),'*_df100',W.file(i)));
    assert(numel(d)==1,'file not found for %s',W.window(i));
    folder = strrep(d.folder,'\','/'); out = fullfile(work,char(W.window(i))); mkdir(out);
    setenv('AG_DEPL',[folder(1:end-1) '*']); setenv('AG_OUT',out); setenv('AG_MASK',char(W.file(i)));
    sg = W.start_segment(i); setenv('AG_SEG0',num2str(sg)); setenv('AG_SEG1',num2str(sg+blk-1));
    evalc('ag_under_test');
    m = dir(fullfile(out,'*','*.mat')); r = load(fullfile(m(1).folder,m(1).name));
    s1 = (sg-1)*N + 1;
    if isempty(r.allSmpPts), det = zeros(0,2); else, det = (r.allSmpPts - s1)/fs; end

    if W.status(i) == "no airguns" || W.status(i) == "unknown", shots = [];
    else, shots = P.time_s(P.window == W.window(i)); end

    % events: normally one per mark; in a hump window, marks closer than join_s form one event
    shots = sort(shots(:)); ev = num2cell(shots');
    hj = Hp.join_s(Hp.window == W.window(i));
    if ~isempty(hj) && ~isempty(shots)
        brk = [true; diff(shots) > hj(1)]; id = cumsum(brk);
        ev = arrayfun(@(q) shots(id == q)',1:id(end),'UniformOutput',false);
    end

    % match each event to the first detection that contains it; a detection that already
    % holds an event marks later events in it as merged. For a hump, a detection anywhere
    % in the hump counts, and every detection in the hump is treated as part of it.
    used = zeros(size(det,1),1); found = 0; merged = 0;
    for q = 1:numel(ev)
        m = ev{q};
        if isscalar(m), k = find(m >= det(:,1)-tol & m <= det(:,2)+tol);
        else, k = find(det(:,2) >= m(1)-tol & det(:,1) <= m(end)+tol); end
        if isempty(k), continue, end
        free = k(used(k)==0);
        if isempty(free), merged = merged + 1; continue, end
        found = found + 1;
        if isscalar(m), used(free(1)) = 1; else, used(k) = 1; end
    end
    res(i,:) = [numel(ev) size(det,1) found numel(ev)-found merged sum(used==0) sum(det(:,2)-det(:,1) >= 5)];

    % measurements of every detection, for looking at what separates real shots from false ones
    if isfield(r,'allTier'), tr = r.allTier; else, tr = nan(size(det,1),1); end
    for k = 1:size(det,1)
        feats(end+1,:) = {W.window(i),det(k,1),det(k,2),used(k),tr(k),r.allCorrVal(k,1), ...
            r.allRmsDet(k)-r.allRmsNBefore(k),r.allRmsDet(k)-r.allRmsNAfter(k), ...
            r.allPpDet(k)-r.allPpNBefore(k),r.allPpDet(k)-r.allPpNAfter(k),r.allDur(k)}; %#ok<AGROW>
    end

    for k = 1:size(det,1), dets(end+1,:) = {W.window(i),det(k,1),det(k,2)}; end %#ok<AGROW>
    if isfield(r,'allCand') && ~isempty(r.allCand)       % every candidate the prototype considered, and what became of it
        tc = r.allCand(:,1) - datenum(r.rawStart(sg,:))*86400;
        for k = 1:numel(tc), cands(end+1,:) = [{W.window(i),tc(k)} num2cell(r.allCand(k,2:end))]; end %#ok<AGROW>
    end

    % shots and detection starts in each of the window's five minutes
    ns = histcounts(shots,0:60:300); nd = histcounts(min(max(det(:,1),0),299.99),0:60:300);
    for k = 1:5, mins(end+1,:) = {W.window(i),k,ns(k),nd(k)}; end %#ok<AGROW>
end
T = table(W.window,W.condition_guess,W.site,W.status,res(:,1),res(:,2),res(:,3),res(:,4),res(:,5),res(:,6),res(:,7), ...
    'VariableNames',{'window','condition','site','status','shots','detections','found','missed','merged','false','long_dets'});
writetable(T,fullfile(setDir,"scores_" + outName + ".csv"));
writetable(cell2table(dets,'VariableNames',{'window','start_s','end_s'}),fullfile(setDir,"detections_" + outName + ".csv"));
writetable(cell2table(feats,'VariableNames',{'window','start_s','end_s','is_shot','tier','score','rms_over_before_dB', ...
    'rms_over_after_dB','pp_over_before_dB','pp_over_after_dB','duration_s'}),fullfile(setDir,"features_" + outName + ".csv"));
if ~isempty(cands)
    writetable(cell2table(cands,'VariableNames',{'window','time_s','score','grade','fail_rms_after','fail_pp_after','fail_rms_before', ...
        'fail_pp_before','fail_duration','duration_s','outcome'}),fullfile(setDir,"candidates_" + outName + ".csv"));
end
writetable(cell2table(mins,'VariableNames',{'window','minute','shots','detections'}),fullfile(setDir,"minutes_" + outName + ".csv"));
if nargout == 0, disp(T); end
a = T.shots > 0;
fprintf('\nAirgun windows (%d): %d shots, %d found (%.0f%%), %d missed (of which %d merged into another detection), %d false detections\n', ...
    sum(a),sum(T.shots(a)),sum(T.found(a)),100*sum(T.found(a))/sum(T.shots(a)),sum(T.missed(a)),sum(T.merged(a)),sum(T.false(a)));
fprintf('No-airgun windows (%d): %d false detections (%.1f per 5-minute window), %d of them 5 s or longer\n', ...
    sum(~a),sum(T.false(~a)),mean(T.false(~a)),sum(T.long_dets(~a)));
end
