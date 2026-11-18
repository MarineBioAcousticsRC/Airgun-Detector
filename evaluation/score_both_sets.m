function score_both_sets
% SCORE_BOTH_SETS  Score the original detector and the shape-based detector
% (settings as written in each file) on test sets 1, 2 and 3, using the
% reviewed marks. Prints one summary line per run.
here = fileparts(mfilename('fullpath')); proto = fullfile(fileparts(here),'xcorr_shape','xcorr_airgun_shape.m');
runs = {'set1', 'reviewed_original', ''
        'set1', 'reviewed_prototype', proto
        'set3', 'reviewed_original', ''
        'set3', 'reviewed_prototype', proto
        'set2', 'reviewed_original', ''
        'set2', 'reviewed_prototype', proto};
for i = 1:size(runs,1)
    T = score_detector(runs{i,2},{},runs{i,3},runs{i,1});
    a = T.shots > 0;
    fprintf('RESULT %-5s %-19s | windows with airguns %d: shots %d, found %d (%.0f%%), merged %d, extra detections %d | windows without airguns %d: false detections %d\n', ...
        runs{i,1},runs{i,2},sum(a),sum(T.shots(a)),sum(T.found(a)),100*sum(T.found(a))/max(1,sum(T.shots(a))),sum(T.merged(a)),sum(T.false(a)),sum(~a),sum(T.false(~a)));
end
end
