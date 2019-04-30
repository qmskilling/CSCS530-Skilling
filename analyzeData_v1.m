clearvars
close all

pathName = '../Data_fastDecay_noSTDP/';

paramDir = dir([pathName,'*Param*.dat']);
paramDir = sortDirByNumber(transpose({paramDir.name}));

connDir = dir([pathName,'*Conn*.dat']);
connDir = transpose({connDir.name});

findTime = ~cellfun('isempty', strfind(connDir,'Times'));
connTime_dir = sortDirByNumber(connDir(findTime));
connDir = sortDirByNumber(connDir(~findTime));

spikeDir = dir([pathName,'*spike*.dat']);
spikeDir = sortDirByNumber(transpose({spikeDir.name}));

weightDir = dir([pathName,'*Weight*.dat']);
weightDir = sortDirByNumber(transpose({weightDir.name}));

paramData = zeros(length(paramDir),4);

outDegreeDistribution = cell(length(paramData),1);
inDegreeDistribution = cell(length(paramData), 1);
largestComponentNodes = cell(length(paramData),1);
outDegreeRatio = zeros(length(paramData),2);
freqData = cell(length(paramData),1);
outDegree_full = cell(length(paramData),1);
inDegree_full = cell(length(paramData),1);
centralityMeasures = cell(length(paramData),1);

for f_ind = 327:length(paramDir)
   
    paramData(f_ind,:) = dlmread([pathName,paramDir{f_ind}]);
    
    spikeData_i = dlmread([pathName,spikeDir{f_ind}]);
    connData_i = dlmread([pathName, connDir{f_ind}]);
    connTimes_i = dlmread([pathName, connTime_dir{f_ind}]);
    weightData_i = dlmread([pathName,weightDir{f_ind}]);
    
    spikeData_i_ro = reorganizeByPosition(spikeData_i);

    [connData_i_ro, connTimes_i_ro] = reorganizeConnsByPosition(connData_i, connTimes_i);
    [adjMat_vs_time, unqTimes] = adjacencyMatrixTrace(connData_i(:,[1 5:end]), connTimes_i(:, [1 5:end]));
    [adjMat_i, adjMat_w_i] = createAdjMat(connData_i(:,[1 5:end]), weightData_i(:,[1 5:end]));
    
    G = digraph(adjMat_i,'omitselfloops');
    G_w = digraph(adjMat_w_i,'omitselfloops');
    
    n = numnodes(G);
    
    uw_BC = centrality(G,'betweenness');
    uw_BC = 2*uw_BC./((n-2)*(n-1));
    inC = centrality(G,'indegree');
    outC = centrality(G,'outdegree');
    
    w_BC = centrality(G_w,'betweenness');
    w_BC = 2*w_BC./((n-2)*(n-1));
    
    centralityMeasures{f_ind} = [inC outC uw_BC w_BC];
    
    a = true;
    
    largestComponentNodes{f_ind} = largestcomponent(adjMat_i);
    fCount = 1;
    
    mpcData = zeros(size(adjMat_vs_time,3)-2,10);
    
    for ind = 2:size(adjMat_vs_time,3)
       
        adjMat_t = adjMat_vs_time(:,:,ind);
        
        outDegree_t = [connData_i(:,1) sum(adjMat_t,2)];
        inDegree_t = [connData_i(:,1) transpose(sum(adjMat_t,1))];
        comp_t = largestcomponent(adjMat_t);
        
        if(length(comp_t) > 1)
           
            nrnsInComp = ismember(1:size(adjMat_i,1), comp_t);
            
            [spikeData_inComp, len_inComp] = sortSpikeData(spikeData_i(:, [1 4:5]),find(nrnsInComp), unqTimes(ind-1:ind));
            %[spikeData_notInComp, len_notInComp] = sortSpikeData(spikeData_i(:,[1 4:5]), find(~nrnsInComp), unqTimes(ind-1:ind));
            
            if(len_inComp > 1)
                mpc_inComp = meanPhaseCoherence(spikeData_inComp);
                %mpc_notInComp = meanPhaseCoherence(spikeData_notInComp);
            end
            
            %figure(1)
            %subplot(6,7,fCount)
            %imagesc(mpc_inComp)
            %c1 = colorbar;
            %colormap('jet')
            %caxis([0 1])
            %c1.Label.String = 'M.P.C.';
            %title(['Growth Step = ',num2str(ind-1)])
            %set(gca,'XTick','')
            %set(gca,'YTick','')
            %subplot(1,2,2)
            %imagesc(mpc_notInComp)
            %c2 = colorbar;
            %colormap('jet')
            %caxis([0 1])
            
            mpc_inComp(1:length(mpc_inComp)+1:end) = [];
            mpc_bins = zeros(length(0.1:0.1:1),1);
            jj = 1;
            
            for ii = 0.1:0.1:1
               
                mpc_bins(jj) = length(find(mpc_inComp >=(ii-0.1) & mpc_inComp <ii))/length(mpc_inComp);
                
                jj = jj + 1;
                
            end
            
            mpcData(fCount,:) = mpc_bins;
            
            fCount =fCount + 1;
            
        end
        
        a = true;
        
    end
    
    imagesc(mpcData)
    c = colorbar;
    colormap('parula')
    caxis([0 1])
    ylabel('Growth Step with Component Size > 1')
    xlabel('M.P.C.')
    c.Label.String = 'P(M.P.C.)';
    set(gca,'XTick',1:10,'XTickLabel',0.1:1)
    
    %imagesc(adjMat_i);
    
    outDegree = [connData_i(:,1) sum(adjMat_i,2)];
    outDegree_full{f_ind} = outDegree;
    inDegree = [connData_i(:,1) transpose(sum(adjMat_i,1))];
    inDegree_full{f_ind} = inDegree;
    %}
   
    freqData_i = calcFreq(spikeData_i(:,[1 4:5]));
    freqData{f_ind} = freqData_i;
    
%{
    outgoing_pairs = get_freqDegreePairs(outDegree, freqData_i);
    incoming_pairs = get_freqDegreePairs(inDegree, freqData_i);
    
    outDegreeDistribution{f_ind} = getDegreeDistribution(outDegree);
    inDegreeDistribution{f_ind} = getDegreeDistribution(inDegree);
    
    useData = connData_i(:,4) > 0;
    
    outDegreeRatio(f_ind,:) = [mean(outDegree(useData,2)./connData_i(useData,4)) std(outDegree(useData,2)./connData_i(useData,4))/sqrt(length(connData_i(useData,1)))];
    %}
    a = true;
    
end

save('graphMeasures_slowDecay.mat','centralityMeasures','freqData');

%save('metaData_4.mat','paramData','outDegreeDistribution','inDegreeDistribution','largestComponentNodes','outDegreeRatio','freqData');
%save('metaData_4-degreesFull.mat','outDegree_full','inDegree_full');

function [sortedDir] = sortDirByNumber(fNames)

findUndsc = strfind(fNames,'_');
findHyp = strfind(fNames,'-');

if(length(find(~cellfun('isempty', findHyp))) == 0)
    
    fileNumber = arrayfun(@(x) str2num(fNames{x}(findUndsc{x}+1:end-4)), 1:length(findUndsc),'UniformOutput',false);
    fileNumber = cell2mat(transpose(fileNumber));
    
else
    
    fileNumber = arrayfun(@(x) str2num(fNames{x}(findUndsc{x}+1:findHyp{x}-1)), 1:length(findUndsc),'UniformOUtput',false);
    fileNumber = cell2mat(transpose(fileNumber));
end

[x,y] = sort(fileNumber);

sortedDir = fNames(y);

end

function [data_ro] = reorganizeByPosition(data)

uniqueID = unique(data(:,1));
findID = arrayfun(@(x) data(:,1)==uniqueID(x), 1:length(uniqueID),'UniformOutput',false);
position = arrayfun(@(x) data(find(findID{x},1,'first'), 2:3), 1:length(uniqueID),'UniformOutput',false);
position = cell2mat(transpose(position));

[x,y] = sortrows(position);
uniqueID_ro = uniqueID(y);
tempData = data;

for u_ind = 1:length(uniqueID_ro)
   
    findID = data(:,1) == uniqueID_ro(u_ind);
    
    tempData(findID,1) = u_ind;
    
end

data_ro = tempData(:,[1 4:5]);

a = true;


end

function [connList, connTimes_list] = reorganizeConnsByPosition(conns, connTimes)

uniqueID = unique(conns(:,1));
findID = arrayfun(@(x) conns(:,1)==uniqueID(x), 1:length(uniqueID),'UniformOutput',false);
position = arrayfun(@(x) conns(findID{x}, 2:3), 1:length(uniqueID),'UniformOutput',false);
position = arrayfun(@(x) position{x}(1,:), 1:length(position),'UniformOutput',false);
position = cell2mat(transpose(position));

[x,y] = sortrows(position);
uniqueID_ro = uniqueID(y);
newIDs = zeros(length(uniqueID_ro),1);

for u_ind = 1:length(uniqueID_ro)
   
    findID = conns(:,1) == uniqueID_ro(u_ind);
    
    newIDs(u_ind) = u_ind;
    
end

IDMap = [uniqueID_ro newIDs];

%Replaces IDs in conn List
connList = conns(:,[1 5:end]);
connTimes_list = connTimes(:, [1 5:end]);

for u_ind = 1:size(connList,1)
   
    connList_i = connList(u_ind,2:end);
    [x,y] = sort(connList_i);
    connList_i = connList_i(y);
    
    connList_i(find(connList_i)) = IDMap(ismember(IDMap(:,1), connList_i(find(connList_i))),2);
    
    connList(u_ind,2:end) = connList_i;
    connList(u_ind,1) = IDMap(ismember(IDMap(:,1), connList(u_ind,1)),2);
    connTimes_list(u_ind,1) = IDMap(ismember(IDMap(:,1), connTimes_list(u_ind,1)),2);
    
end

[x,y] = sort(connList(:,1));

connList = connList(y,:);
connTimes_list = connTimes_list(y,:);

a = true;

end

function [adjMat, adjMat_w] = createAdjMat(connData, wghtData)

adjMat = zeros(size(connData,1));
adjMat_w = zeros(size(connData,1));

for ind = 1:size(connData,1)
   
    adjMat(connData(ind,1),connData(ind,connData(ind,2:end)>0)) = 1;
    adjMat_w(connData(ind,1), connData(ind, connData(ind,2:end)>0)) = wghtData(ind,connData(ind,2:end)>0);
    
end

end

function [freq] = calcFreq(spkData)

[x,y] = sort(spkData(:,2));
sData = spkData(y,:);

uniqueID = unique(spkData(:,1));
findID = arrayfun(@(x) spkData(:,1)==uniqueID(x), 1:length(uniqueID),'UniformOutput',false);
ISI = arrayfun(@(x) diff(spkData(findID{x},3)), 1:length(uniqueID),'UniformOutput',false);
freq = arrayfun(@(x) [spkData(find(findID{x},1,'first'),1:2) mean(1000./ISI{x}) std(1000./ISI{x})/sqrt(length(ISI{x}))], 1:length(uniqueID),'UniformOutput',false);
freq = cell2mat(transpose(freq));

[x,y] = sort(freq(:,1));
freq = freq(y,:);

a = true;

end

function [freqDegree_pairs] = get_freqDegreePairs(degree, freq)


    freqDegree_pairs = zeros(size(degree,1),2);
    
    for ind = 1:size(degree,1)
       
        findID = degree(ind,1) == freq(:,1);
        
        if(isempty(find(findID,1)))
           continue; 
        end
        
        freqDegree_pairs(ind,:) = [freq(find(findID,1),3) degree(ind,2)];
        
        a = true;
        
    end
    
end

function [distribution] = getDegreeDistribution(degree)

uniqueDegree = unique(degree(:,2));
histDegree = histc(degree(:,2), uniqueDegree)./length(degree);

distribution = [uniqueDegree histDegree];

end

function [phaseData] = meanPhaseCoherence(spkData)

uniqueID = unique(spkData(:,1));

phaseData = zeros(length(uniqueID));

for  u_ind = 1:length(uniqueID)
   
    data_i = spkData(spkData(:,1)==uniqueID(u_ind),3);
    
    for u_jnd = u_ind+1:length(uniqueID)
        
        data_j = spkData(spkData(:,1)==uniqueID(u_jnd),3);
        
        phases_ij = getPhases(data_i, data_j);
        %phases_ij = [mean(phases_ij) std(phases_ij)/sqrt(length(phases_ij))];
        phases_ji = getPhases(data_j, data_i);
        %phases_ji = [mean(phases_ji) std(phases_ji)/sqrt(length(phases_ji))];
        
        phaseData(u_ind, u_jnd) = mean(phases_ij);
        phaseData(u_jnd, u_ind) = mean(phases_ji);
        
        a = true;
        
    end
    
    a = true;
    
end

%phaseData = phaseData.*8*atan(1.0);
%phaseData(1:length(phaseData)+1:end) = [];

%meanPhase = mean(phaseData)/(8*atan(1.));
%meanPhase = mean(phaseData);

a = true;

    function [phases] = getPhases(d_i, d_j)
       
        phases = [];
        
        for ii = 1:length(d_i)-1
           
            findData = d_j >= d_i(ii) & d_j < d_i(ii+1);
            
            if(isempty(find(findData,1)))
                phases = [phases; nan];
                continue;
            end
            
            phases = [phases; (d_j(findData)-d_i(ii))./(d_i(ii+1)-d_i(ii))];
            
        end
        
        phases(isnan(phases)) = [];
        a = true;
        
    end

a = true;

end

function [adjMat_t, uniqueTimes] = adjacencyMatrixTrace(connList, connTimes)

uniqueTimes = unique(connTimes(:,2:end));

adjMat_t = zeros(size(connList,1), size(connList,1), length(uniqueTimes));

for t_ind = 2:length(uniqueTimes)
   
    if(uniqueTimes(t_ind) == 0)
        continue;
    end
    
    adjMat_t(:,:,t_ind) = adjMat_t(:,:,t_ind-1);
    
    for ind = 1:size(connList,1)
        
        adjMat_t(connList(ind,1), connList(ind,connTimes(ind,2:end)==uniqueTimes(t_ind)), t_ind) = 1;
        
    end

    
end


a = true;
function [B] = largestcomponent(A)
%%%NOT ORIGINAL
%%% SOURCE: https://www.mathworks.com/matlabcentral/fileexchange/30926-largest-component

n=length(A);
for i=1:n
    mz{i}=find(A(i,:));
end
x(1:n)=0;
z=0;
k=0;
for i=1:n
    if x(i)==0;
        z=z+1;
        clear v
        v(1)=i;
        while nnz(v)>0
                    x(v(1))=z;
                    k=union(k,v(1));
                    b=setdiff(mz{v(1)},k);
                    v=setdiff(v,v(1));
                    v=union(v,b);         
        end
    end
end
c(1:max(x))=0;
for i=1:max(x)
    c(i)=length(find(x==i)); 
end
cm=find(c==max(c));
    
B=find(x==cm(1));
end

end

function [B] = largestcomponent(A)
%%%NOT ORIGINAL
%%% SOURCE: https://www.mathworks.com/matlabcentral/fileexchange/30926-largest-component

n=length(A);
for i=1:n
    mz{i}=find(A(i,:));
end
x(1:n)=0;
z=0;
k=0;
for i=1:n
    if x(i)==0;
        z=z+1;
        clear v
        v(1)=i;
        while nnz(v)>0
                    x(v(1))=z;
                    k=union(k,v(1));
                    b=setdiff(mz{v(1)},k);
                    v=setdiff(v,v(1));
                    v=union(v,b);         
        end
    end
end
c(1:max(x))=0;
for i=1:max(x)
    c(i)=length(find(x==i)); 
end
cm=find(c==max(c));
    
B=find(x==cm(1));
end

function [connRatio] = getConnRatio(connData)

expectedConns = connData(:,4);

test = arrayfun(@(x) find(connData(x,5:end)), 1:size(connData,1),'UniformOutput',false);
numConns = 0;

a = true;


end

function [spkData_selected, lenIDs] = sortSpikeData(spkData, nrns, times)

uniqueIDs = unique(spkData(:,1));
uniqueIDs = uniqueIDs(ismember(uniqueIDs, nrns));
spkData_IDs = arrayfun(@(x) spkData(spkData(:,1)==uniqueIDs(x),:), 1:length(uniqueIDs),'UniformOutput',false);

spkData_selected = cell2mat(transpose(spkData_IDs));
spkData_selected = spkData_selected(spkData_selected(:,3)>=times(1)*0.1 & spkData_selected(:,3) < times(2)*0.1,:);

lenIDs = length(uniqueIDs);


a = true;

end