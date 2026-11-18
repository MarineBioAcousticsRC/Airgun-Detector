function plot_methods(nameA,nameB)
% PLOT_METHODS  Figures that show how the original detector and the
% prototype each look at the same stretch of audio.
%   plot_methods('baseline','prototype_frozen')
% Writes figures/M_template.png and one three-panel figure per example:
%   (a) the low-passed waveform and its envelope, with the reviewed shots
%   (b) the original detector's match curve and threshold
%   (c) the prototype's shape score and its two thresholds
% Detections come from detections_<name>.csv written by score_detector.

fs = 2000; N = 150000; base = fileparts(mfilename('fullpath')); here = fullfile(base,'set1'); root = '//frosty.ucsd.edu/GofMX_Decimated_8';
outDir = fullfile(base,'figures'); if ~exist(outDir,'dir'), mkdir(outDir); end
W = readtable(fullfile(here,'windows.csv'),'Delimiter',',','TextType','string');
P = readtable(fullfile(here,'truth.csv'),'Delimiter',',','TextType','string');
A = readtable(fullfile(here,"detections_" + nameA + ".csv"),'Delimiter',',','TextType','string');
B = readtable(fullfile(here,"detections_" + nameB + ".csv"),'Delimiter',',','TextType','string');
colA = [0.84 0.37 0.00]; colB = [0.00 0.45 0.70]; ink = [0.15 0.15 0.15]; grey = [0.62 0.65 0.68];

% settings of the two methods (as in xcorr_airgun_df100.m and xcorr_airgun_shape.m)
c2_offset = 0.15; minCorrThresh = 7e15; shapeSure = 0.65; shapeWeak = 0.40; envSmooth_s = 0.1;
tp = load(fullfile(fileparts(base),'xcorr_original','air_template_df100.mat')); template = tp.DATA(:)'; lt = numel(template);
envT = abs(hilbert(template));
ds = round(fs/50); nSm = round(envSmooth_s*fs);
et = movmean(envT,nSm); et = et(1:ds:end); M = numel(et); tz = et-mean(et); tNorm = sqrt(sum(tz.^2));
[~,iPulse] = max(et); pulseOff_s = (iPulse-1)*ds/fs;
[Bf,Af] = butter(5,150/(fs/2),'low');

% ---- the template ----
fig = figure('Visible','off','Position',[50 50 1500 420],'Color','w'); ax = axes(fig,'Position',[0.06 0.16 0.91 0.72]); hold(ax,'on');
tt = (0:lt-1)/fs; sc = max(abs(template));
plot(ax,tt,template/sc,'Color',grey,'LineWidth',0.6);
plot(ax,tt,movmean(envT,nSm)/sc,'Color',ink,'LineWidth',1.8);
legend(ax,{'template waveform','its envelope (what both methods compare against)'},'Location','northeast','Box','off');
xlabel(ax,'Time (s)'); ylabel(ax,'Amplitude (relative)'); xlim(ax,[0 lt/fs]); set(ax,'TickDir','out','Box','off');
title(ax,'The 6-second airgun template (air\_template\_df100.mat)','FontWeight','normal');
exportgraphics(fig,fullfile(outDir,'M_template.png'),'Resolution',120); close(fig);

%          window  segment-in-window  file name      title
ex = {'W07', 1, 'M_fast',  'Fast shooting (W07, CE\_03): shots 6.5 s apart'
      'W10', 2, 'M_ship',  'Ship noise, no airguns (W10, AC\_02)'
      'W22', 1, 'M_faint', 'Faint shots (W22, Y2B\_01)'};
for x = 1:size(ex,1)
    w = W(W.window==ex{x,1},:); k = ex{x,2};
    d = dir(fullfile(root,"GOM_" + w.site,'*_df100',w.file)); f = fullfile(d.folder,d.name);
    s1 = (w.start_segment-1)*N + 1;                     % first sample of the window
    a = s1 + (k-1)*N;                                    % first sample of the segment shown
    t0 = (k-1)*75; tlim = [t0 t0+75];
    % the segment alone, as the original detector reads it
    y = double(audioread(f,[a a+N-1],'native')); y = y-mean(y);
    yF = filtfilt(Bf,Af,y); env = abs(hilbert(yF.'));
    [c,lg] = xcorr(env,envT); c = c(lg>=0); c2 = c.^2;
    med = median(c2(1:end-lt)); thr = max(med*(1+c2_offset),minCorrThresh);
    tc = t0 + (0:numel(c2)-1)/fs + lt/2/fs;             % shifted half a template, as the detector does
    % the segment with context, as the prototype reads it
    pre = lt; post = 2*lt; yy = double(audioread(f,[a-pre a+N-1+post],'native')); yy = yy-mean(yy);
    envS = movmean(abs(hilbert(filtfilt(Bf,Af,yy))),nSm); e = envS(1:ds:end)';
    num = conv(e,fliplr(tz),'valid'); q1 = conv(e,ones(1,M),'valid'); q2 = conv(e.^2,ones(1,M),'valid');
    r = num./(tNorm*sqrt(max(q2-q1.^2/M,eps)));
    tr = t0 - pre/fs + (0:numel(r)-1)*ds/fs + pulseOff_s; % time of the pulse each score refers to

    shots = P.time_s(P.window==ex{x,1}); if contains(w.condition_guess,'no airguns'), shots = []; end
    da = A(A.window==ex{x,1},:); db = B(B.window==ex{x,1},:);

    fig = figure('Visible','off','Position',[50 50 1500 900],'Color','w');
    tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
    title(tl,ex{x,4},'FontSize',12);

    ax1 = nexttile(tl); hold(ax1,'on'); te = t0 + (0:N-1)/fs;
    plot(ax1,te,yF,'Color',grey,'LineWidth',0.4); plot(ax1,te,movmean(env,nSm),'Color',ink,'LineWidth',1.3);
    yl = ylim(ax1); h = plot(ax1,shots,yl(2)*0.93*ones(size(shots)),'v','MarkerSize',7,'MarkerFaceColor',ink,'MarkerEdgeColor','none');
    if isempty(shots), legend(ax1,{'low-passed waveform','envelope'},'Location','northeastoutside','Box','off');
    else, legend(ax1,{'low-passed waveform','envelope','reviewed shots'},'Location','northeastoutside','Box','off'); end
    ylabel(ax1,'Amplitude (counts)'); title(ax1,'(a) The audio, low-passed at 150 Hz','FontWeight','normal','FontSize',10);

    ax2 = nexttile(tl); hold(ax2,'on'); set(ax2,'YScale','log');
    lo = max(min(c2(c2>0))*0.8,thr/1e4); hi = max(max(c2),thr)*2.5; hd = [];
    for j = 1:height(da), hd = patch(ax2,[da.start_s(j) da.end_s(j) da.end_s(j) da.start_s(j)],[lo lo hi hi],colA,'FaceAlpha',0.22,'EdgeColor','none'); end
    h1 = plot(ax2,tc,max(c2,lo),'Color',ink,'LineWidth',1.1); h2 = yline(ax2,thr,'-','Color',colA,'LineWidth',1.6);
    ylim(ax2,[lo hi]); ylabel(ax2,'Match (squared correlation)');
    if thr == minCorrThresh, tn = 'threshold (the fixed floor, 7e15)'; else, tn = 'threshold (15% above the segment median)'; end
    if isempty(hd), legend(ax2,[h1 h2],{'match curve',tn},'Location','northeastoutside','Box','off');
    else, legend(ax2,[h1 h2 hd],{'match curve',tn,'detections'},'Location','northeastoutside','Box','off'); end
    title(ax2,'(b) Original detector: a detection is any stretch where the match curve is above the threshold','FontWeight','normal','FontSize',10);

    ax3 = nexttile(tl); hold(ax3,'on'); hd = [];
    for j = 1:height(db), hd = patch(ax3,[db.start_s(j) db.end_s(j) db.end_s(j) db.start_s(j)],[-0.6 -0.6 1 1],colB,'FaceAlpha',0.25,'EdgeColor','none'); end
    h1 = plot(ax3,tr,r,'Color',ink,'LineWidth',1.1);
    h2 = yline(ax3,shapeSure,'-','Color',colB,'LineWidth',1.6); h3 = yline(ax3,shapeWeak,'--','Color',colB,'LineWidth',1.3);
    ylim(ax3,[-0.6 1]); ylabel(ax3,'Shape score (-1 to 1)'); xlabel(ax3,'Time in window (s)');
    names = {'shape score',sprintf('%.2f: accepted on shape alone',shapeSure),sprintf('%.2f: accepted only if in a train',shapeWeak)};
    if isempty(hd), legend(ax3,[h1 h2 h3],names,'Location','northeastoutside','Box','off');
    else, legend(ax3,[h1 h2 h3 hd],[names {'detections'}],'Location','northeastoutside','Box','off'); end
    title(ax3,'(c) Prototype: a candidate shot is each peak of the shape score','FontWeight','normal','FontSize',10);

    set([ax1 ax2 ax3],'XLim',tlim,'XTick',t0:5:t0+75,'TickDir','out','Box','off'); linkaxes([ax1 ax2 ax3],'x');
    exportgraphics(fig,fullfile(outDir,[ex{x,3} '.png']),'Resolution',120); close(fig);
    fprintf('%s: threshold %.3g (median %.3g, floor in use %d); max match %.3g; shape score max %.2f\n',ex{x,1},thr,med,thr==minCorrThresh,max(c2),max(r));
end
end
