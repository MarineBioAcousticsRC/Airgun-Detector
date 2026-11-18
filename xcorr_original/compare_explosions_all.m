% get a list of airgun detector output files.
fList = dir('E:\Data\Explosions\*.mat'); % these files generated with airgun_df100_saveSnips.m

maxLag = 2000; % maximum lag steps for xcorr
myDist = [];
Ncomp = 200;% number of subsequent detections to compare to.
for iFile = 1%:nFiles
    load(fullfile(fList(iFile).folder, fList(iFile).name))
    nDets = size(expAll,1);
    myDist = zeros(nDets);

    for i1 = 1:nDets
        %Step through detections, and compare each detection to next N
        %detections
        myDet1 = expAll(i1,:);
        for i2 = (i1+1): (min(i1+Ncomp,nDets))
            myDet2 = expAll(i2,:);
            
            [r,lags] = xcorr(myDet1,myDet2,maxLag);
            
            [maxxcorr,maxXcorrIdx] = max(r);
            
            if false
                myX1 = 1:length(myDet1);
                myX2 = [1:length(myDet2)]+lags(maxXcorrIdx);
                myEndIdx = min([myX1(end),myX2(end),20000]);
                figure(11);clf
                subplot(2,1,1)
                plot(myDet1);hold on;plot(myDet2)
                subplot(2,1,2)
                plot(myDet1);hold on;plot(myX2,myDet2)
            end
            
            
                myDist(i1,i2) = maxxcorr;
           
            
        end
        fprintf('Done with detection %0.0f of %0.0f\n',i1, nDets)
    end
end

%%


mySim = myDist/10^7;
minSimilarity = prctile(mySim(mySim>0),75,'all');
mySim(mySim<=minSimilarity)=0;
% for iS = 1:size(mySim,1)
%     [~,mIdx]= maxk(mySim(iS,:),1);
%     setZero = setdiff(1:size(mySim,2),mIdx);
%     mySim(iS,setZero) = 0;
% end
p.maxCWiterations = 20;
p.plotFlag = 1;
allRowIdx = 1:size(mySim,1);
clusterIDtemp = ct_run_CW_cluster(allRowIdx,mySim,p);
clustBins = 0:max(clusterIDtemp);
counts = histc(clusterIDtemp,clustBins);
keepClust = find(counts >= 5);
clustNums = clustBins(keepClust);
nodeSet = {};
clusterID = nan(size(clusterIDtemp))';
iC = 1;
if ~isempty(keepClust)
    for i4 = 1:length(keepClust)
        nodeIndex = clusterIDtemp==clustNums(i4);
        nodeSet{i4,1} = allRowIdx(nodeIndex);
        clusterID(nodeIndex) = iC;
        iC = iC+1;
    end
end
clf
imagesc(MSN(isnan(clusterID),:))
imagesc(MSN(~isnan(clusterID),:))

MTT = btKeepAll(:,4);
MSN = expAll;
MPP = 20*log10(max(MSN,[],2)+abs(min(MSN,[],2)));
TPWSname = strrep(fList(iFile).name,'.mat','_TPWS1.mat');
save(TPWSname,'MSN','MTT','MPP','-v7.3')

zID(:,1) = MTT(:,1);
zID(:,2) = clusterID;
zID(isnan(clusterID),:) = [];

IDname = strrep(fList(iFile).name,'.mat','_ID1.mat');
save(IDname,'zID','-v7.3')

% figure(20);clf
% G = graph(mySim(1:400,1:400),'upper','omitselfloops');
% h = plot(G,'layout','force','nodelabel',[1:400]);
% h.EdgeColor = [.5,.5,.5];
% mySimOrig = mySim;
% zID(:,1) = timeStack(:,1);
% zID(:,2) = clusterID; % ID zero will break detEdit, must delete these before finishing up.
% checkedIdx = zeros(size(timeStack(:,1))); % this will track which detections have
% % been evaluated as part of a train.
% mySeq1 = 1;
% checkedIdx(mySeq1,1) = 1;
% IDnum = 1;
% while mySeq1< size(mySim, 1)
%     mySeqAll = mySeq1;
%     
%     iCount = 2;
%     while max(mySeqAll) <size(mySim, 1)
%         nextOne = find(mySim(mySeq1,:)>0)+max(mySeqAll);
%         
%         if isempty(nextOne)&& size(mySeqAll,1)>1
%             break
%         elseif max(mySeqAll)==size(mySim, 1)
%             break
%         elseif isempty(nextOne)
%             break
%         end
%         mySeqAll(iCount,1) = nextOne;
%         mySeq1 = nextOne;
%         iCount = iCount+1;
%         checkedIdx(mySeq1,1) = 1;
%                 
%     end
%     mySeq1 = find(checkedIdx==0,1,'first');
%     checkedIdx(mySeq1,1) = 1;
%     if size(mySeqAll,1)>1
%         zID(mySeqAll,2) = IDnum;
%         if IDnum==15
%             %reset ID num if it hits max
%             IDnum = 1;
%         else
%             IDnum = IDnum+1;
%         end
%     end
%     mySeqAll = mySeq1;
% 
% end
% clf(2)
% for iS = 1:length(mySeqAll)
%     plot((double(audioStack{mySeqAll(iS)})./max(double(audioStack{mySeqAll(iS)})))+...
%         ((timeStack(iS,1)-(timeStack(1,1)))*6*60*24),'k')
%     hold on
% end


