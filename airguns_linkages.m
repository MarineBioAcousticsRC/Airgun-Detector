% relationlist

% MMATCH: each row represents an explosion, columns contain indices of
% similar subsequent explosions. 20 explosions are reviewed, the best 4 are
% given, in order of time.
% MTT: Each row represents an explosion, start and end times of explosion are given
% MSIM: Similarity score of explosions in MMATCH
% 
MMATCH = [];
MTT = [];
MSIM = [];

for iR = 1:size(mySim, 1)
    [bestMatchVal,bestMatchIdx] = maxk(mySim(iR,iR+(1:10)),4);
    [~,iX] = sort(bestMatchIdx,'ascend');
    MMATCH(iR,:) = iR+bestMatchIdx(iX);
    MTT(iR,:) = timeStack(iR,:);
    MSIM(iR,:) = bestMatchVal(iX);
    CSN{iR,1} = audioStack{};
end

% grab the series of detections most similar to the starting event (mySeq1)
% this could go sideways if there's a bad match, so we probably want to skip 
% explosions that are below some similarity score.
mySeq1 = 1;
iCount = 2;
while mySeq1 <size(mySim, 1)
    mySeq1(iCount) = min(MMATCH(mySeq1(iCount-1),:));
    
    iCount = iCount+1;
end

figure(11);clf
for iExp = 1:length(mySeq1)
    plot()
    
end