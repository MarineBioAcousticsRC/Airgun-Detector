function plot_comparison(nameA,labelA,nameB,labelB,setName,nameC,labelC)
% PLOT_COMPARISON  Spectrogram of each test window with the reviewed shots
% and the detections from two or three scored runs drawn above it.
%   plot_comparison('reviewed_original','original','reviewed_prototype','shape-based','set2')
%   plot_comparison('reviewed_original','original','reviewed_prototype','shape-based','set2','network','network')
% reads detections_<name>.csv from the set's folder (written by score_detector,
% or copied there from the network's output).

if nargin < 5 || isempty(setName), setName = 'set1'; end
three = nargin >= 7;
fs = 2000; N = 150000; blk = 4; here = fileparts(fullfile(fileparts(mfilename('fullpath')),setName,'x')); root = '//frosty.ucsd.edu/GofMX_Decimated_8';
W = readtable(fullfile(here,'windows.csv'),'Delimiter',',','TextType','string');
P = readtable(fullfile(here,'truth.csv'),'Delimiter',',','TextType','string');
if exist(fullfile(here,'review.csv'),'file')
    Rv = readtable(fullfile(here,'review.csv'),'Delimiter',',','TextType','string');
    [~,ix] = ismember(W.window,Rv.window); W.condition_guess = Rv.status(ix); W = W(W.condition_guess ~= "dropped",:);
end
A = readtable(fullfile(here,"detections_" + nameA + ".csv"),'Delimiter',',','TextType','string');
B = readtable(fullfile(here,"detections_" + nameB + ".csv"),'Delimiter',',','TextType','string');
if three, C = readtable(fullfile(here,"detections_" + nameC + ".csv"),'Delimiter',',','TextType','string'); if ~isstring(C.window), C.window = string(C.window); end, end
outDir = fullfile(here,'plots_compare'); if ~exist(outDir,'dir'), mkdir(outDir); end
colA = [0.84 0.37 0.00]; colB = [0.00 0.45 0.70]; colC = [0.00 0.62 0.45]; top = 3 + three;
for i = 1:height(W)
    d = dir(fullfile(root,"GOM_" + W.site(i),'*_df100',W.file(i))); sg = W.start_segment(i);
    y = double(audioread(fullfile(d.folder,d.name),[(sg-1)*N+1,(sg+blk-1)*N],'native')); y = y-mean(y);
    [~,F,T,Pw] = spectrogram(y,hann(2000),1600,2000,fs); keep = F<=250; PdB = 10*log10(Pw(keep,:));
    noAir = contains(W.condition_guess(i),'no airguns') || W.condition_guess(i) == "unknown";
    if noAir, shots = []; else, shots = P.time_s(P.window==W.window(i)); end
    a = A(A.window==W.window(i),:); b = B(B.window==W.window(i),:);

    fig = figure('Visible','off','Position',[50 50 1700 720],'Color','w');
    % top strip: reviewed shots and the two sets of detections
    ax1 = axes(fig,'Position',[0.085 0.77-0.03*three 0.86 0.15+0.03*three]); hold(ax1,'on');
    plot(ax1,shots,top*ones(size(shots)),'v','MarkerSize',7,'MarkerFaceColor',[0.15 0.15 0.15],'MarkerEdgeColor','none');
    for k = 1:height(a), plot(ax1,[a.start_s(k) a.end_s(k)],[top-1 top-1],'-','Color',colA,'LineWidth',5); plot(ax1,a.start_s(k),top-1,'|','Color',[0.4 0.15 0],'MarkerSize',11,'LineWidth',1.2); end
    for k = 1:height(b), plot(ax1,[b.start_s(k) b.end_s(k)],[top-2 top-2],'-','Color',colB,'LineWidth',5); plot(ax1,b.start_s(k),top-2,'|','Color',[0 0.2 0.35],'MarkerSize',11,'LineWidth',1.2); end
    labs = {sprintf('%s (%d)',labelB,height(b)),sprintf('%s (%d)',labelA,height(a))};
    if three   % the network gives a start time only, so it is drawn as a tick without a bar
        c = C(C.window==W.window(i),:);
        plot(ax1,c.start_s,ones(height(c),1),'|','Color',colC,'MarkerSize',13,'LineWidth',2.2);
        labs = [{sprintf('%s (%d)',labelC,height(c))},labs];
    end
    if W.condition_guess(i) == "unknown", lab1 = 'reviewed shots (not sure)'; elseif noAir, lab1 = 'reviewed shots (none)'; else, lab1 = sprintf('reviewed shots (%d)',numel(shots)); end
    set(ax1,'YLim',[0.4 top+0.6],'YTick',1:top,'YTickLabel',[labs,{lab1}], ...
        'XLim',[0 300],'XTick',0:10:300,'XTickLabel',[],'TickDir','out','Box','off','YGrid','on','GridColor',[0.85 0.85 0.85],'GridAlpha',1,'FontSize',9);
    title(ax1,sprintf('%s   |   %s   |   %s   |   starts %s        (bars run from detection start to end; tick marks the start)', ...
        W.window(i),W.condition_guess(i),strrep(W.site(i),'_','\_'),W.start_time(i)),'FontWeight','normal','FontSize',10);
    % spectrogram
    ax2 = axes(fig,'Position',[0.085 0.09 0.86 0.66-0.03*three]);
    imagesc(ax2,T,F(keep),PdB); axis(ax2,'xy'); colormap(ax2,flipud(gray(256))); clim(ax2,[prctile(PdB(:),5) prctile(PdB(:),99.9)]);
    set(ax2,'XLim',[0 300],'XTick',0:10:300,'TickDir','out','YLim',[0 250]); xlabel(ax2,'Time (s)'); ylabel(ax2,'Frequency (Hz)');
    for bnd = 75:75:225, xline(ax2,bnd,':','Color',[0.3 0.3 0.3]); end
    exportgraphics(fig,fullfile(outDir,W.window(i) + ".png"),'Resolution',110); close(fig);
end
for s = 1:ceil(height(W)/8)
    ids = (s-1)*8 + (1:8); ids = ids(ids<=height(W));
    ims = arrayfun(@(k) imresize(imread(fullfile(outDir,W.window(k) + ".png")),0.5),ids,'UniformOutput',false);
    imwrite(imtile(ims,'GridSize',[4 2],'BorderSize',4,'BackgroundColor','w'),fullfile(outDir,sprintf('contact_sheet_%d.png',s)));
end
end
