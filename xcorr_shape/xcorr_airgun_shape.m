function  xcorr_airgun_shape
% XCORR_AIRGUN_SHAPE  Shape-based two-stage airgun detector.
% Run it on a deployment with run_xcorr_shape.
%
% Stage 1 - find candidate shots by shape.
%   The low-passed signal envelope is compared with the template envelope
%   using the correlation coefficient (-1 to 1), which measures how alike
%   the two shapes are regardless of how loud the signal is. Each peak of
%   that curve above parm.shapeThresh is one candidate shot. Peaks must be
%   at least parm.minSep_s apart, so neighbouring shots stay separate.
%
%   Two additions to stage 1:
%   - A pulse that stands far above the background in the plain (level-
%     dependent) match curve is accepted whatever its shape. This keeps
%     loud shots whose shape differs from the template, such as several
%     arrivals overlapping (parm.loudProm).
%   - Peaks between parm.shapeWeak and parm.shapeThresh are "weak"
%     candidates, unless they stand at least parm.clearPP_dB above the
%     noise before them. A weak candidate is kept only if it is part of a train:
%     other detections must sit at a steady spacing before and/or after it
%     (parm.reg*). Airgun shots repeat at a fixed interval; noise does not.
%
% Stage 2 - mark each candidate's start and end, join the arrivals of one
%   shot into a single event, then apply the same signal-versus-noise tests
%   as xcorr_airgun_df100 (rms and peak-to-peak against the noise before and
%   after) to the event. An event is rejected when parm.minFlags or more of
%   the tests fail.
%   Joining: a shot often reaches the recorder as two or three arrivals a
%   few seconds apart. Candidates no more than parm.groupGap_s apart, with
%   no deep dip in the envelope between them (parm.groupValley_dB), become
%   one event that runs from the first arrival to the end of the last. The
%   noise before the event is measured ahead of its first arrival.
%
% Differences from xcorr_airgun_df100 besides stage 1:
%   - No fixed floor (minCorrThresh) and no median-based threshold.
%   - Each segment is read with a few seconds of the neighbouring segments
%     when they are continuous in time, so shots near a segment boundary
%     are seen whole. A shot is reported by the segment its pulse falls in.
%   - Start and end are where the smoothed envelope falls parm.edge_dB
%     below the pulse's peak.
%   - Duration is not used to reject short events: with start and end
%     measured this way a sharp pulse is naturally well under a second.
%     Only events of parm.durLong_s or longer count as a failed test.
%   - allCorrVal holds [shape score at the peak, median shape score of the
%     segment] instead of squared-correlation values.
% Output variables and file layout are otherwise the same.

global PARAMS
%parm stores parameters, to be used in TethysXML output
parm.shapeThresh = 0.65; %correlation coefficient at or above this: accepted on shape alone
parm.shapeWeak = 0.4; %between shapeWeak and shapeThresh: kept only if part of a regular train (set equal to shapeThresh to switch off)
parm.clearPP_dB = 8; %a weak candidate standing this far (peak-to-peak, dB) above the noise before it is accepted without a train (Inf switches off)
parm.loudProm = 1; %a match-curve peak standing this many times the segment median above its surroundings is accepted whatever its shape (Inf switches off)
parm.regMin_s = 4; %shortest shot spacing considered when looking for a train (s)
parm.regMax_s = 40; %longest shot spacing considered (s)
parm.regTol_s = 0.6; %how far a shot may sit from its expected place in a train (s)
parm.regSupport = 3; %number of other detections that must line up with a weak candidate
parm.minSep_s = 2.5; %minimum time between candidate shots (s)
parm.envSmooth_s = 0.1; %smoothing applied to envelopes before comparing shapes (s)
parm.edge_dB = 10; %start/end of a shot: where the smoothed envelope is this far below its peak
parm.nSamples = 1000; %number of noise samples to be pulled out
parm.nBuff = 2000; %number of samples to buffer between end of detection and noise sample.
parm.rmsAS = 1; %rms noise after signal <rmsAS (dB) difference counts as a failed test
parm.rmsBS = 1; %rms noise before signal
parm.ppAS = 1; %pp noise after signal <ppAS (dB) difference counts as a failed test
parm.ppBS = 1; %pp noise before signal
parm.durLong_s = 15; %duration >= durLong_s (s) counts as a failed test
parm.durShort_s = 0; %duration <= durShort_s (s) counts as a failed test (0 switches this off)
parm.maxTail_s = 5; %an arrival's end is looked for no later than this long after its peak (s)
parm.groupGap_s = 3.5; %arrivals no more than this far apart may be joined into one event (0 switches joining off)
parm.groupValley_dB = Inf; %...but only if the envelope between them dips less than this below the smaller arrival (Inf: ignore dips)
parm.minFlags = 2; %number of failed tests needed to reject a candidate
overWritePrevious = 1;
deplList = dir('\\frosty.ucsd.edu\GofMX_Decimated_8\GOM_*');
outDisk  = 'N:\tempAirguns';
templateFile = fullfile(fileparts(fileparts(mfilename('fullpath'))),'xcorr_original','air_template_df100.mat');
SearchFileMask = {'*df100.x.wav';'*d100.x.wav'};

templateStruct = load(templateFile);
template = templateStruct.DATA(:)';
len_template = length(template);

for iDepl = 1:length(deplList)

    BaseDir = fullfile(deplList(iDepl).folder,deplList(iDepl).name);
    outputDir = fullfile(outDisk,[deplList(iDepl).name,'_AirgunDetections']);
    if ~exist(outputDir,'dir'), mkdir(outputDir); end
    [~, FileList, PathList] = utFindFiles(SearchFileMask, {BaseDir}, 1);

    for fidx = 1:size(FileList,1)
        file = FileList{fidx};
        fend = strfind(file,'.df100');
        if isempty(fend)
            fend = strfind(file,'.x.wav');
        end
        newFile = fullfile(outputDir,[file(1:fend-1),'.mat']);
        if ~overWritePrevious && exist(newFile,'file')
            disp('Output file already exists, overwrite false, skipping.')
            continue
        end
        filepath = fullfile(PathList{fidx},file);
        display(['calculating ',file,'; file ',num2str(fidx),'/',num2str(size(FileList,1))])
        [rawStart,rawDur,fs] = readxwavhd(filepath);
        inc = size(rawStart,1);

        % where each segment sits in the file, in samples
        bps = PARAMS.samp.byte;
        segFirst = (PARAMS.xhd.byte_loc - PARAMS.xhd.byte_loc(1))/bps + 1;
        segLen = PARAMS.xhd.byte_length/bps;
        segTime = datenum(rawStart);
        % contiguous(k) is true when segment k+1 follows segment k without a gap in time
        contiguous = abs(diff(segTime(:))*86400 - segLen(1:end-1)'/fs) < 0.01;
        fileAudioInfo = audioinfo(filepath);
        fileSampleNum = fileAudioInfo.TotalSamples;

        % low-pass filter, as in xcorr_airgun_df100
        [B,A] = butter(5, 150/(fs/2),'low');

        % template envelope, smoothed and reduced to about 50 points per second
        ds = max(1,round(fs/50));
        nSmooth = max(1,round(parm.envSmooth_s*fs));
        et = movmean(abs(hilbert(template)),nSmooth); et = et(1:ds:end);
        M = numel(et); tz = et - mean(et); tNorm = sqrt(sum(tz.^2));
        [~,iPulse] = max(et); pulseOff = (iPulse-1)*ds;   % samples from template start to its pulse

        allSmpPts = []; allExp = []; allCorrVal = []; allDur = [];
        allRmsNBefore = []; allRmsNAfter = []; allRmsDet = [];
        allPpNBefore = []; allPpNAfter = []; allPpDet = [];
        allTier = []; allPulse = []; allCand = []; candRow = [];

        for idx = 1:inc
            a = segFirst(idx); n = segLen(idx);
            if a + n - 1 > fileSampleNum
                warning('File has fewer than expected samples. Skipping this segment.');
                continue
            end
            % extra audio from the neighbouring segments, when continuous
            pre = 0; post = 0;
            if idx > 1 && contiguous(idx-1), pre = min(len_template,segLen(idx-1)); end
            if idx < inc && contiguous(idx), post = min(2*len_template,segLen(idx+1)); end
            last = min(a + n - 1 + post, fileSampleNum);

            y = double(audioread(filepath,[a-pre last],'native'));
            y = y - mean(y);
            yFilt = filtfilt(B,A,y);
            envS = movmean(abs(hilbert(yFilt)),nSmooth);
            envG = movmean(envS,round(0.5*fs));   % smoother copy, for judging dips between arrivals
            e = envS(1:ds:end)';
            if numel(e) < M + 3, continue, end

            % ---- stage 1: shape score and its peaks ----
            num = conv(e,fliplr(tz),'valid');
            s1 = conv(e,ones(1,M),'valid'); s2 = conv(e.^2,ones(1,M),'valid');
            r = num./(tNorm*sqrt(max(s2 - s1.^2/M,eps)));
            minDist = min(round(parm.minSep_s*fs/ds),numel(r)-2);
            pk = []; loc = [];
            if any(r >= parm.shapeWeak)
                [pk,loc] = findpeaks(r,'MinPeakHeight',parm.shapeWeak,'MinPeakDistance',minDist);
            end
            tier = 1 + (pk < parm.shapeThresh);     % 1 = accepted outright, 2 = weak (needs a train)
            if isfinite(parm.loudProm)
                % loud pulses of any shape, from the plain match curve
                c2 = conv(e,fliplr(et),'valid').^2;
                [~,locL,~,prL] = findpeaks(c2,'MinPeakDistance',minDist);
                for q = locL(prL >= parm.loudProm*median(c2))
                    [dmin,k] = min(abs(loc - q));
                    if isempty(loc) || dmin > minDist
                        loc(end+1) = q; pk(end+1) = r(q); tier(end+1) = 1; %#ok<AGROW>
                    else
                        tier(k) = 1;
                    end
                end
                [loc,order] = sort(loc); pk = pk(order); tier = tier(order);
            end
            if isempty(loc), continue, end
            p = (loc-1)*ds + 1 + pulseOff;          % estimated pulse position in y
            medScore = median(r);

            % ---- stage 2a: start, end and peak of each candidate arrival ----
            nc = numel(p); cs = zeros(1,nc); ce = zeros(1,nc); cpk = zeros(1,nc);
            for j = 1:nc
                lo = max(1,p(j)-fs); hi = min(numel(y),p(j)+2*fs);
                [pv,ii] = max(envS(lo:hi)); ipk = lo + ii - 1;
                edge = pv*10^(-parm.edge_dB/20);
                k = find(envS(lo:ipk) < edge,1,'last');
                if isempty(k), s = lo; else, s = lo + k; end
                limit = min(numel(y),ipk + round(parm.maxTail_s*fs));
                if j < nc, limit = max(ipk+1,min(limit,p(j+1)-fs)); end
                k = find(envS(ipk:limit) < edge,1,'first');
                if isempty(k), en = limit; else, en = ipk + k - 2; end
                cs(j) = s; ce(j) = max(en,s+1); cpk(j) = ipk;
            end

            % ---- stage 2b: join the arrivals of one shot into a single event ----
            grp = ones(1,nc);
            for j = 2:nc
                i1 = min(cpk(j-1),cpk(j)); i2 = max(cpk(j-1),cpk(j));
                dip_dB = 20*log10(min(envG(i1),envG(i2))/max(min(envG(i1:i2)),eps));
                joined = (i2 - i1) <= parm.groupGap_s*fs && dip_dB < parm.groupValley_dB;
                grp(j) = grp(j-1) + ~joined;
            end

            % ---- stage 2c: signal-versus-noise tests on each event ----
            for g = 1:grp(end)
                mem = find(grp == g);
                ipk = cpk(mem(1));                                % the first arrival places the event in time
                if ipk <= pre || ipk > pre + n, continue, end     % belongs to a neighbouring segment
                s = min(cs(mem)); en = max(ce(mem));
                score = max(pk(mem)); tierG = min(tier(mem));

                yDet = yFilt(s:en);
                yNAfter = yFilt(min(en+parm.nBuff,numel(yFilt)+1):min(en+parm.nBuff+parm.nSamples,numel(yFilt)));
                yNBefore = yFilt(max(s-parm.nSamples-parm.nBuff,1):(s-1-parm.nBuff));

                ppDetSeg = 20*log10(max(yDet)-min(yDet));
                rmsDetSeg = 20*log10(sqrt(mean(yDet.^2)));
                if isempty(yNAfter), ppNAfterSeg = NaN; rmsNAfterSeg = NaN;
                else, ppNAfterSeg = 20*log10(max(yNAfter)-min(yNAfter)); rmsNAfterSeg = 20*log10(sqrt(mean(yNAfter.^2))); end
                if isempty(yNBefore), ppNBeforeSeg = NaN; rmsNBeforeSeg = NaN;
                else, ppNBeforeSeg = 20*log10(max(yNBefore)-min(yNBefore)); rmsNBeforeSeg = 20*log10(sqrt(mean(yNBefore.^2))); end
                durSeg = (en - s)/fs;

                % the tests of xcorr_airgun_df100 (missing noise windows never fail a test)
                flags = [abs(rmsDetSeg - rmsNAfterSeg) < parm.rmsAS, ...
                         abs(ppDetSeg - ppNAfterSeg) < parm.ppAS, ...
                         (rmsDetSeg - rmsNBeforeSeg) < parm.rmsBS, ...
                         (ppDetSeg - ppNBeforeSeg) < parm.ppBS, ...
                         durSeg >= parm.durLong_s || durSeg <= parm.durShort_s];
                nFail = sum(flags);
                % a record of every candidate, for diagnosis: pulse time (s), shape score, grade,
                % the five tests (1 = failed), duration, and the outcome
                % (1 = rejected by the tests, 2 = rejected by the train rule, 3 = kept)
                allCand = [allCand; segTime(idx)*86400 + (ipk-pre)/fs, score, tierG, flags, durSeg, 1 + 2*(nFail < parm.minFlags)];
                if nFail >= parm.minFlags, continue, end
                if tierG == 2 && (ppDetSeg - ppNBeforeSeg) >= parm.clearPP_dB, tierG = 1; end

                % positions relative to this segment, then to the file and to time
                rel = [s en] - pre;
                allSmpPts = [allSmpPts; rel + a - 1];
                allExp = [allExp; segTime(idx) + (rel/fs)/(60*60*24)];
                allCorrVal = [allCorrVal; score medScore];
                allDur = [allDur; durSeg];
                allRmsNBefore = [allRmsNBefore; rmsNBeforeSeg];
                allRmsNAfter = [allRmsNAfter; rmsNAfterSeg];
                allRmsDet = [allRmsDet; rmsDetSeg];
                allPpNBefore = [allPpNBefore; ppNBeforeSeg];
                allPpNAfter = [allPpNAfter; ppNAfterSeg];
                allPpDet = [allPpDet; ppDetSeg];
                allTier = [allTier; tierG];
                allPulse = [allPulse; segTime(idx)*86400 + (ipk-pre)/fs];
                candRow = [candRow; size(allCand,1)];
            end
        end

        % weak candidates stay only if they are part of a regularly spaced train
        if ~isempty(allExp)
            keep = true(size(allTier));
            for i = find(allTier == 2)'
                keep(i) = fits_train(allPulse,i,parm);
            end
            allSmpPts = allSmpPts(keep,:); allExp = allExp(keep,:); allCorrVal = allCorrVal(keep,:);
            allDur = allDur(keep); allRmsNBefore = allRmsNBefore(keep); allRmsNAfter = allRmsNAfter(keep);
            allRmsDet = allRmsDet(keep); allPpNBefore = allPpNBefore(keep); allPpNAfter = allPpNAfter(keep);
            allPpDet = allPpDet(keep); allTier = allTier(keep);
            allCand(candRow(~keep),end) = 2;
        end

        bt = [];
        if ~isempty(allExp)
            bt(:,1) = allSmpPts(:,1)-5000;
            bt(:,2) = allSmpPts(:,2)+5000;
            bt(:,3) = zeros(size(allSmpPts,1),1);
            bt(:,4) = allExp(:,1);
            bt(:,5) = allExp(:,2);
        end
        save(newFile,'allSmpPts','allExp','allCorrVal','allDur',...
            'allRmsNBefore','allRmsNAfter','allRmsDet','allPpNBefore',...
            'allPpNAfter','allPpDet','allTier','allCand','rawStart','rawDur','parm','bt','-v7.3');
    end
end
end

function ok = fits_train(t,i,parm)
% True if detection i sits in a regularly spaced train: for some spacing T
% between parm.regMin_s and parm.regMax_s, at least parm.regSupport other
% detections lie at t(i) + k*T for k = -3..3. t is in seconds.
d = t - t(i); ok = false;
d = d(abs(d) <= 3*parm.regMax_s + 2*parm.regTol_s);   % only nearby detections can matter; keeps whole-file runs fast
for j = find(abs(d) >= parm.regMin_s & abs(d) <= parm.regMax_s)'
    T = abs(d(j)); n = 0;
    for k = [-3 -2 -1 1 2 3]
        n = n + any(abs(d - k*T) <= parm.regTol_s*(1 + 0.25*(abs(k)-1)));
    end
    if n >= parm.regSupport, ok = true; return, end
end
end
