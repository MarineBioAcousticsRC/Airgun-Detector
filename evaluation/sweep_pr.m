function sweep_pr
% SWEEP_PR  Precision and recall of the original detector and the prototype
% over a range of settings, from strict to lenient, on set 1 and set 2.
%
%   recall    = shots found / shots marked
%   precision = detections that hold a shot / all detections
%               (detections in windows without airguns all count against it)
%
% Original detector: the fixed floor on the squared correlation
% (minCorrThresh) is swept; everything else stays as written. With the
% floor at 0 only the "15% above the segment median" rule is left, so the
% last points lower that offset as well.
%
% Prototype: it has no single threshold, so its four acceptance settings
% are moved together, from strict to lenient (shapeWeak, shapeThresh,
% loudProm, clearPP_dB). The row marked current is the file as written.
%
% Results go to pr_curve.csv; plot_pr.py draws the curves together with the
% network's.

here = fileparts(mfilename('fullpath')); proto = fullfile(fileparts(here),'xcorr_shape','xcorr_airgun_shape.m');
%        floor     offset   current
origS = {'1e17',   '0.15',  0
         '5e16',   '0.15',  0
         '2e16',   '0.15',  0
         '7e15',   '0.15',  1
         '3e15',   '0.15',  0
         '1e15',   '0.15',  0
         '3e14',   '0.15',  0
         '1e14',   '0.15',  0
         '0',      '0.15',  0
         '0',      '0.05',  0
         '0',      '0',     0};
%         weak    thresh   loudProm  clearPP  current
protoS = {'0.8',  '0.95',  'Inf',    'Inf',   0
          '0.7',  '0.9',   'Inf',    'Inf',   0
          '0.6',  '0.85',  '4',      '16',    0
          '0.5',  '0.75',  '2',      '12',    0
          '0.4',  '0.65',  '1',      '8',     1
          '0.3',  '0.55',  '0.75',   '6',     0
          '0.2',  '0.45',  '0.5',    '4',     0
          '0.1',  '0.35',  '0.25',   '2',     0};
sets = {'set1','set2'}; rows = {};
for s = 1:2
    for i = 1:size(origS,1)
        T = score_detector('pr_tmp',{'minCorrThresh',origS{i,1}; 'parm.c2_offset',origS{i,2}},'',sets{s});
        rows(end+1,:) = add(T,'original',sprintf('floor %s offset %s',origS{i,1:2}),i,origS{i,3},s); %#ok<AGROW>
    end
    for i = 1:size(protoS,1)
        T = score_detector('pr_tmp',{'parm.shapeWeak',protoS{i,1}; 'parm.shapeThresh',protoS{i,2}; ...
            'parm.loudProm',protoS{i,3}; 'parm.clearPP_dB',protoS{i,4}},proto,sets{s});
        rows(end+1,:) = add(T,'prototype',sprintf('weak %s shape %s loud %s clear %s',protoS{i,1:4}),i,protoS{i,5},s); %#ok<AGROW>
    end
end
R = cell2table(rows,'VariableNames',{'detector','setting','step','current','set','shots','detections','found','spare'});
writetable(R,fullfile(here,'pr_curve.csv'));
for s = 1:2                                  % remove the temporary score files
    delete(fullfile(here,sets{s},'*_pr_tmp.csv'));
end
end

function r = add(T,det,label,step,cur,s)
r = {det,label,step,cur,sprintf('set%d',s),sum(T.shots),sum(T.detections),sum(T.found),sum(T.false)};
fprintf('PR %-9s set%d %-42s found %d of %d, detections %d, without a shot %d\n',det,s,label,r{8},r{6},r{7},r{9});
end
