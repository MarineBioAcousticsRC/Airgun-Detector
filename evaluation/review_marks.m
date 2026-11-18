function review_marks(setName,selftest)
% REVIEW_MARKS  Point-and-click editor for the shot marks in a test set.
%
%   review_marks          first test set  (this folder)
%   review_marks('set2')  second test set (the set2 folder)
%
% Shows each five-minute window as a spectrogram (top) and the 10-120 Hz
% envelope (bottom) with the shot marks drawn on both. Changes are saved
% automatically to truth.csv (the marks) and review.csv (what each window
% is) in the set's folder.
%
%   Mouse   left click          add a mark at that time
%           right click         remove the nearest mark (within 3 s)
%   Keys    right / left arrow  next / previous window
%           k                   marks are right: window has airguns
%           x                   window has NO airguns (marks are ignored)
%           d                   drop the window from the test set (e.g. instrument noise)
%           n                   not sure whether these are airguns (scored as no airguns, but flagged)
%           m                   airguns are there, but the shots cannot be marked one by one
%                               (kept as an airgun window, left out of the shot-by-shot scores)
%           u                   back to "not reviewed"
%           c                   clear all marks in this window
%           r                   reset this window to the automatic proposals
%           z                   zoom in to 75 s around the pointer;  a  shows all 300 s
%           h                   show / hide hints
%           q                   save and close
%
% Hints: if the set has a hints.csv (see make_hints.py), each hint is drawn as a
% blue triangle along the bottom of the spectrogram. A hint marks a place where a
% detector fired and there is no mark. It is only a prompt to look there.
%
% The second argument is for automated checking only.

if nargin < 1 || isempty(setName), setName = 'set1'; end
if nargin < 2, selftest = false; end
fs = 2000; N = 150000; blk = 4; root = '//frosty.ucsd.edu/GofMX_Decimated_8';
here = fileparts(mfilename('fullpath')); setDir = fullfile(here,setName);
truthFile = fullfile(setDir,'truth.csv'); reviewFile = fullfile(setDir,'review.csv');

W = readtable(fullfile(setDir,'windows.csv'),'Delimiter',',','TextType','string');
Pk = readtable(fullfile(setDir,'picks.csv'),'Delimiter',',','TextType','string');
nW = height(W);

% marks: start from truth.csv if it exists, otherwise from the automatic proposals
marks = cell(nW,1); proposals = cell(nW,1);
for i = 1:nW, proposals{i} = sort(Pk.time_s(Pk.window == W.window(i)))'; end
if exist(truthFile,'file')
    Tr = readtable(truthFile,'Delimiter',',','TextType','string');
    for i = 1:nW, marks{i} = sort(Tr.time_s(Tr.window == W.window(i)))'; end
else
    marks = proposals;
end
% status of each window
status = repmat("not reviewed",nW,1);
if exist(reviewFile,'file')
    Rv = readtable(reviewFile,'Delimiter',',','TextType','string');
    [tf,ix] = ismember(W.window,Rv.window); status(tf) = Rv.status(ix(tf));
elseif exist(truthFile,'file')          % a set reviewed before this tool existed
    status(:) = "airguns"; status(contains(W.condition_guess,'no airguns')) = "no airguns";
end

% optional hints: places where a detector fired with no mark nearby
hints = cell(nW,1); showHints = true;
if exist(fullfile(setDir,'hints.csv'),'file')
    Hn = readtable(fullfile(setDir,'hints.csv'),'Delimiter',',','TextType','string');
    if ~isstring(Hn.window), Hn.window = string(Hn.window); end
    for i = 1:nW, hints{i} = Hn.time_s(Hn.window == W.window(i))'; end
end

cache = cell(nW,1); cur = 1; pointerTime = 150; hMarks = gobjects(0);
[bb,ab] = butter(4,[10 120]/(fs/2),'bandpass');

if selftest, vis = 'off'; else, vis = 'on'; end
fig = figure('Name','Review shot marks','NumberTitle','off','Color','w','Visible',vis, ...
    'Units','normalized','Position',[0.04 0.08 0.92 0.82],'MenuBar','none','ToolBar','none', ...
    'WindowButtonDownFcn',@onClick,'WindowButtonMotionFcn',@onMove,'KeyPressFcn',@onKey,'CloseRequestFcn',@onClose);
ax1 = axes(fig,'Position',[0.05 0.36 0.93 0.56]);
ax2 = axes(fig,'Position',[0.05 0.11 0.93 0.20]);
uicontrol(fig,'Style','text','Units','normalized','Position',[0.02 0.005 0.96 0.045],'BackgroundColor','w', ...
    'HorizontalAlignment','left','FontSize',9,'String', ...
    ['Left click: add mark    Right click: remove nearest mark    Arrows: next / previous window    ' ...
     'k: marks are right    x: no airguns    n: not sure    m: airguns but shots unclear    d: drop window    u: not reviewed    c: clear marks    ' ...
     'r: reset to proposals    z: zoom at pointer    a: show all    h: hints on/off    q: save and close']);
showWindow(1);

if selftest
    % scripted actions, used to check the tool without a person
    n0 = numel(marks{1}); addMark(123.4); assert(numel(marks{1}) == n0+1);
    removeNear(123.0); assert(numel(marks{1}) == n0);
    setStatus("no airguns"); showWindow(2); setStatus("airguns"); addMark(10);
    pointerTime = 200; zoomAt(); xl = xlim(ax1); assert(abs(diff(xl)-75) < 1e-6); xlim([ax1 ax2],[0 300]);
    showWindow(2); saveAll();
    if ischar(selftest) || isstring(selftest), exportgraphics(fig,selftest,'Resolution',80); end
    delete(fig);
    fprintf('selftest ok: %d windows, files written to %s\n',nW,setDir);
end

    function loadWindow(i)
        if ~isempty(cache{i}), return, end
        title(ax1,sprintf('Loading %s ...',W.window(i)),'FontWeight','normal'); drawnow;
        d = dir(fullfile(root,"GOM_" + W.site(i),'*_df100',W.file(i)));
        sg = W.start_segment(i);
        y = double(audioread(fullfile(d(1).folder,d(1).name),[(sg-1)*N+1,(sg+blk-1)*N],'native')); y = y-mean(y);
        [~,F,T,P] = spectrogram(y,hann(2000),1600,2000,fs); keep = F <= 250;
        c.T = T; c.F = F(keep); c.PdB = 10*log10(P(keep,:));
        e = movmean(abs(hilbert(filtfilt(bb,ab,y))),0.25*fs); e = e(1:fs/20:end);
        c.et = (0:numel(e)-1)/20; c.edB = 20*log10(e+eps); c.edB = c.edB - movmedian(c.edB,30*20);
        cache{i} = c;
    end

    function showWindow(i)
        cur = min(max(i,1),nW); loadWindow(cur); c = cache{cur};
        cla(ax1); imagesc(ax1,c.T,c.F,c.PdB); axis(ax1,'xy'); colormap(ax1,flipud(gray(256)));
        clim(ax1,[prctile(c.PdB(:),5) prctile(c.PdB(:),99.9)]); hold(ax1,'on');
        set(ax1,'XLim',[0 300],'YLim',[0 250],'XTick',0:10:300,'TickDir','out'); ylabel(ax1,'Frequency (Hz)');
        cla(ax2); plot(ax2,c.et,c.edB,'Color',[0.15 0.15 0.15]); hold(ax2,'on');
        yline(ax2,4,':','Color',[0.5 0.5 0.5]);
        set(ax2,'XLim',[0 300],'XTick',0:10:300,'TickDir','out','YLim',[min(-6,min(c.edB)) max(12,max(c.edB)+2)]);
        xlabel(ax2,'Time (s)'); ylabel(ax2,{'10-120 Hz envelope','(dB above background)'});
        hMarks = gobjects(0); drawMarks();
    end

    function drawMarks()
        delete(hMarks(isgraphics(hMarks))); hMarks = gobjects(0);
        m = marks{cur}; col = [0.84 0.37 0];
        if ismember(status(cur),["no airguns","dropped","unknown","airguns unclear"]), col = [0.7 0.7 0.7]; end
        y2 = ylim(ax2);
        for k = 1:numel(m)
            hMarks(end+1) = plot(ax1,[m(k) m(k)],[0 232],'-','Color',[col 0.45],'LineWidth',1); %#ok<AGROW>
            hMarks(end+1) = plot(ax1,m(k),238,'v','MarkerSize',8,'MarkerFaceColor',col,'MarkerEdgeColor','w'); %#ok<AGROW>
            hMarks(end+1) = text(ax1,m(k),246,num2str(k),'FontSize',7,'HorizontalAlignment','center','Clipping','on'); %#ok<AGROW>
            hMarks(end+1) = plot(ax2,[m(k) m(k)],y2,'-','Color',[col 0.6],'LineWidth',1); %#ok<AGROW>
        end
        if showHints
            for hh = hints{cur}
                hMarks(end+1) = plot(ax1,hh,7,'^','MarkerSize',8,'MarkerFaceColor',[0 0.45 0.70],'MarkerEdgeColor','w'); %#ok<AGROW>
            end
        end
        switch status(cur)
            case "airguns",      st = sprintf('AIRGUNS, %d shots',numel(m));
            case "no airguns",   st = 'NO AIRGUNS (marks ignored)';
            case "dropped",      st = 'DROPPED from the test set';
            case "unknown",      st = 'NOT SURE (scored as no airguns, flagged)';
            case "airguns unclear", st = 'AIRGUNS, SHOTS UNCLEAR (left out of shot-by-shot scores)';
            otherwise,           st = sprintf('not reviewed yet, %d proposed shots',numel(m));
        end
        title(ax1,sprintf('%s  (%d of %d)   |   %s   |   starts %s   |   %s        [reviewed %d of %d]', ...
            W.window(cur),cur,nW,strrep(W.site(cur),'_','\_'),W.start_time(cur),st,sum(status ~= "not reviewed"),nW), ...
            'FontWeight','normal');
    end

    function addMark(t)
        marks{cur} = sort([marks{cur} round(t,2)]);
        if status(cur) ~= "airguns", status(cur) = "airguns"; end
        drawMarks(); saveAll();
    end

    function removeNear(t)
        m = marks{cur}; if isempty(m), return, end
        [dmin,k] = min(abs(m - t)); if dmin > 3, return, end
        m(k) = []; marks{cur} = m;
        if status(cur) == "not reviewed", status(cur) = "airguns"; end
        drawMarks(); saveAll();
    end

    function setStatus(s)
        status(cur) = s; drawMarks(); saveAll();
    end

    function confirmSeen()
        % leaving a window without changing it counts as accepting it as shown
        if status(cur) == "not reviewed"
            if isempty(marks{cur}), status(cur) = "no airguns"; else, status(cur) = "airguns"; end
            saveAll();
        end
    end

    function zoomAt()
        lo = min(max(pointerTime - 37.5,0),225); xlim([ax1 ax2],[lo lo+75]);
        set([ax1 ax2],'XTick',0:5:300);
    end

    function t = pointerIn()
        % time under the pointer if it is inside either plot, otherwise empty
        t = [];
        for ax = [ax1 ax2]
            cp = get(ax,'CurrentPoint'); xl = xlim(ax); yl = ylim(ax);
            if cp(1,1) >= xl(1) && cp(1,1) <= xl(2) && cp(1,2) >= yl(1) && cp(1,2) <= yl(2), t = cp(1,1); return, end
        end
    end

    function onMove(~,~)
        t = pointerIn(); if ~isempty(t), pointerTime = t; end
    end

    function onClick(~,~)
        t = pointerIn(); if isempty(t), return, end
        if strcmp(get(fig,'SelectionType'),'normal'), addMark(t); else, removeNear(t); end
    end

    function onKey(~,evt)
        switch evt.Key
            case 'rightarrow', confirmSeen(); showWindow(cur+1);
            case 'leftarrow',  confirmSeen(); showWindow(cur-1);
            case 'k', setStatus("airguns");
            case 'x', setStatus("no airguns");
            case 'd', setStatus("dropped");
            case 'n', setStatus("unknown");
            case 'm', setStatus("airguns unclear");
            case 'u', setStatus("not reviewed");
            case 'c', marks{cur} = []; drawMarks(); saveAll();
            case 'r', marks{cur} = proposals{cur}; drawMarks(); saveAll();
            case 'z', zoomAt();
            case 'a', xlim([ax1 ax2],[0 300]); set([ax1 ax2],'XTick',0:10:300);
            case 'h', showHints = ~showHints; drawMarks();
            case 'q', onClose();
        end
    end

    function saveAll()
        % truth.csv holds the marks of every window that has airguns
        ft = fopen(truthFile,'w'); fprintf(ft,'window,time_s,source\n');
        for w = 1:nW
            if ismember(status(w),["no airguns","dropped","unknown","airguns unclear"]), continue, end
            for tm = marks{w}, fprintf(ft,'%s,%.2f,reviewed\n',W.window(w),tm); end
        end
        fclose(ft);
        fr = fopen(reviewFile,'w'); fprintf(fr,'window,status,shots\n');
        for w = 1:nW
            n = numel(marks{w}); if ismember(status(w),["no airguns","dropped","unknown","airguns unclear"]), n = 0; end
            fprintf(fr,'%s,%s,%d\n',W.window(w),status(w),n);
        end
        fclose(fr);
    end

    function onClose(~,~)
        saveAll(); delete(fig);
    end
end
