% This code compares airgun detections based on waveform crosscorrelation
% (pairwise distance is also possible but commented out as of 3/31/22)

% It compares each detection with the next N detections, and picks the
% most similar detection, it then uses a clustering algorithm to associate
% those events into linked trains. Trains are usually between 3 and 100
% events long and end when the similarity gets weak. There is no rule
% preventing trains from merging.

% The code produces zID files in which the same id number (1-15) is applied
% to a series of sequential supposedly-related detections. Once the series
% breaks down, the ID label will end. The idea is that these could be
% displayed in detEdit as different colors.

% INPUTS: Output files from airgun_df100_saveSnips.m

% REQUIREMENTS: Triton clustering remora needs to be on your path.

% OUTPUTS:
% - TPWS1 file: NOTE this only includes MTT, MSN, MPP, MSP, f. Also
% includes startIdx and endIdx which tell you where in the original
% waveform the MSN signal was extracted. True time of the MSN should be
% derivable from MTT(1)+ (startIdx*sample rate*datenum conversion)
% - ID1 file: zID and mySpID (should control colors in detEdit).

% Known issues:
% - Hard-coded indices are used to try to ignore the padding
% on the waveforms, padding wasn't helpful for waveform comparison. These
% hard coded values might cause issues if detections are not roughly
% centered in the snip or are very short.

% - Times are of the original snips in the airgun detector output. MTT
% actually has 2 columns, start and end of the snip. Airgun detections are
% of variable lengths, but the MSN representations have been truncated, so
% the time has changed. I think the exact times of MSN can be determined
% by the startIdx and endIdx variables.

% - 


%% Comparison step
% get a list of airgun detector output files.
fList = dir('P:\AirgunClustering\MC10\*.mat'); % these files generated with airgun_df100_saveSnips.m

maxLag = 8000; % maximum lag steps for xcorr
mySim = [];
Ncomp = 10;% number of subsequent detections to compare to.
nfft = 2000;
overlap = 50;
noverlap = round((overlap/100)*nfft);
fs = 2000;
maxColors = 15;
nFiles = length(fList);
for iFile = 1:nFiles
    load(fullfile(fList(iFile).folder, fList(iFile).name))
    nDets = length(audioStack);
    mySim = zeros(nDets);
    MSP = [];
    MPP = [];
    for i1 = 1:nDets
        %Step through detections, and compare each detection to next N
        %detections
        myDet1 = audioStack{i1};
        window = hanning(length(myDet1));
        [MSP(i1,:),f] = pwelch(myDet1,window,noverlap,nfft,fs);
        MPP(i1,:) = 20*log10(max(myDet1)-(min(myDet1)));
        
        for i2 = (i1+1): (min(i1+Ncomp,nDets))
            myDet2 = audioStack{i2};
            
            [r,lags] = xcorr(myDet1(7001:end-7000),myDet2(7001:end-7000),maxLag);
            
            [maxxcorr,maxXcorrIdx] = max(r);
            
            if false % plotting, turned off, but might be useful to check alignments
                myX1 = 1:length(myDet1);
                myX2 = [1:length(myDet2)]+lags(maxXcorrIdx);
                myEndIdx = min([myX1(end),myX2(end),20000]);
                clf(1)
                subplot(3,1,1)
                plot(myDet1);hold on;plot(myDet2)
                title('Original waveforms')
                subplot(3,1,2)
                plot(myDet1);hold on;plot(myX2,myDet2)
                title('Aligned waveforms')
                subplot(3,1,3)
                title('Truncated & aligned waveforms') % this only matters if you're doing pairwise distances.
                plot(myDet1(myX1>=maxLag&myX1<=myEndIdx))
                hold on;plot(myDet2(myX2>=maxLag&myX2<=myEndIdx))
            end
            
            try
                mySim(i1,i2) = maxxcorr;
                % Alternative distance metric:
                % mySim(i1,i2) = -exp(pdist2(myDet1(myX1>=8000&myX1<=myEndIdx)',myDet2(myX2>=8000&myX2<=myEndIdx)','euclidean'));
                % using this might require some adjustments to the
                % similarity percentile thing.
            catch
                fprintf('comparison failed for %d,%d\n',i1,i2)
            end
            
        end
        fprintf('Done with detection %0.0f of %0.0f\n',i1, nDets)
    end
    
    
    %% Clustering step:
    
    
    % Prune similarity matrix
    mySim = mySim/10^8;
    minSimilarity = prctile(mySim(mySim>0),50,'all');
    mySim(mySim<=minSimilarity)=0;
    for iS = 1:size(mySim,1)
        [~,mIdx]= maxk(mySim(iS,:),1);
        setZero = setdiff(1:size(mySim,2),mIdx);
        mySim(iS,setZero) = 0;
    end
    disp('Done computing similarity')

    p.maxCWiterations = 15; % making this smaller will result in smaller clusters.
    p.plotFlag = 1; % plots will only show if there are fewer than 3000 detections, the network plotting code is slow
    
    allRowIdx = 1:size(mySim,1);
    % Run clustering, this requires the cluster remora in Triton to be on the
    % path
    clusterIDtemp = ct_run_CW_cluster(allRowIdx,mySim,p);
    disp('Done clustering')

    % renumber clusters so that they will show up as colors in detEdit.
    clustBins = 0:max(clusterIDtemp);
    counts = histc(clusterIDtemp,clustBins);
    keepClust = find(counts >= 3);
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
    figure(12);clf
    
    plot(timeStack(:,1),mod(clusterID,maxColors)+1,'.')
    xlabel('Time')
    ylabel('ID number')

    
    MTT = timeStack;
    MSN = zeros(size(audioStack,1),16001);
    startIdx = [];
    endIdx = [];
    for iSN = 1:length(audioStack)
        [~,maxIdx] = max(audioStack{iSN});
        startIdx(iSN,1) = max(maxIdx-8000,1);
        endIdx(iSN,1) = min(maxIdx+8000,length(audioStack{iSN}));
        embedLen = length(startIdx(iSN,1):endIdx(iSN,1));
        embedIdx = (1:embedLen)+floor(abs((embedLen-16001)/2));
        MSN(iSN,embedIdx) = audioStack{iSN}(startIdx(iSN,1):endIdx(iSN,1))';
    end
    TPWSname = strrep(fList(iFile).name,'.mat','_TPWS1.mat');
    
    % save output file
    save(TPWSname,'MSN','MTT','MSP','MSN','f','startIdx','endIdx','-v7.3')
    IDname = strrep(fList(iFile).name,'.mat','_ID1.mat');
    mySpID = cellstr(num2str([1:maxColors]'));
    save(IDname,'zID','mySpID','-v7.3')
end

% figure(20);clf
% G = graph(mySim(1:400,1:400),'upper','omitselfloops');
% h = plot(G,'layout','force','nodelabel',[1:400]);
% h.EdgeColor = [.5,.5,.5];
% mySimOrig = mySim;
% zID(:,1) = timeStack(:,1);
% zID(:,2) = 0; % ID zero will break detEdit, must delete these before finishing up.
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


