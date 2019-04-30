clearvars
close all

load('metaData_3.mat');
load('graphMeasures_fastDecay.mat');
%load('metaData_2-degreesFull.mat');

%centralityMeasures = [(indegree Centrality) (outdegree centrality) 
%                       (unweighted betweenness centrality) (weighted betweenness centrality)];

uniqueRadius = unique(paramData(:,1));

for r_ind = 5:length(uniqueRadius)
   
    findR = paramData(:,1)==uniqueRadius(r_ind);
    
    pData_r = paramData(findR,3:end);
    centMeas_r = centralityMeasures(findR);
    freqData_r = freqData(findR);
    
    unqOG = unique(pData_r(:,1));
    
    for og_ind = 1:length(unqOG)
       
        findOG = pData_r(:,1) == unqOG(og_ind);
        
        pData_og = pData_r(findOG,2);
        centMeas_og = centMeas_r(findOG);
        freqData_og = freqData_r(findOG);
        
        arrayfun(@(x) plotData(centMeas_og{x}, freqData_og{x}, unqOG(og_ind), pData_og(x),x+10*(og_ind-1)), 1:length(centMeas_og),'UniformOutput',false);
        
        a = true;
        
    end
    
    a = true;
    
end

function [] = plotData(cm, fData, out_HM,in_HM, panel)

IDs = 1:500;
findIDs = ismember(IDs,fData(:,1));

figure(1)
subplot(10,10, panel)
plot(cm(findIDs,2), cm(findIDs,3),'o')
%}

figure(2)
subplot(10,10,panel)
plot(fData(:,3), cm(findIDs,3),'o')

if(panel <= 10)
    figure(1)
    subplot(10,10, panel)
    title(['I. H.M. = ',num2str(in_HM)])
    
    figure(2)
    subplot(10,10, panel)
    title(['I. H.M. = ',num2str(in_HM)])
end

if(panel == 1 || (mod(panel-1,10) == 0 && panel ~= 91))
    figure(1)
    subplot(10,10, panel)
    ylabel({['O. H.M. = ',num2str(out_HM)],''})
    
    figure(2)
    subplot(10,10, panel)
    ylabel({['O. H.M. = ',num2str(out_HM)],''})
elseif(panel == 91)
    figure(1)
    subplot(10,10, panel)    
    ylabel({['O. H.M. = ',num2str(out_HM)],'Betweenness Centrality'})
    
    figure(2)
    subplot(10,10, panel)    
    ylabel({['O. H.M. = ',num2str(out_HM)],'Betweenness Centrality'})
end

if(panel >= 91)
    figure(1)
    subplot(10,10, panel)
    xlabel('Outdegree Centrality')
    
    figure(2)
    subplot(10,10, panel)
    xlabel('Freq [Hz]')
end

a = true;

end