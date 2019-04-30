clearvars
close all

load('metaData_3.mat');
load('metaData_3-degreesFull.mat');

uniqueRadi = unique(paramData(:,1));

for r_ind = length(uniqueRadi)-1:length(uniqueRadi)
   
    findRadi = paramData(:,1) == uniqueRadi(r_ind);
    
    freqData_r = freqData(findRadi);
    inFull_r = inDegree_full(findRadi);
    inDegDist_r = inDegreeDistribution(findRadi);
    comp_r = largestComponentNodes(findRadi);
    outFull_r = outDegree_full(findRadi);
    outDegDist_r = outDegreeDistribution(findRadi);
    degRatio_r = outDegreeRatio(findRadi,:);
    pData_r = paramData(findRadi,3:4);
    
    unique_ihm = unique(pData_r(:,2));
    unique_ohm = unique(pData_r(:,1));
    
    compLength = cellfun('length', comp_r);
    compLength = transpose(reshape(compLength, [length(unique(pData_r(:,2))), length(unique(pData_r(:,1)))]));
    figure(1)
    %subplot(2,3,r_ind)
    imagesc(compLength./500)
    colormap('gray')
    c = colorbar;
    caxis([0 1])
    title(['Conn. Radius =',num2str(uniqueRadi(r_ind))])
    set(gca,'XTick','')
    set(gca,'YTick','')
    
    if(r_ind == 4)
        set(gca,'XTick',1:2:9,'XTickLabel',unique_ihm(1:2:9))
        set(gca,'YTick',1:2:9,'YTickLabel',unique_ohm(1:2:9))
        xlabel('Incoming Half-Maximum')
        ylabel('Outgoing Half-Maximum')
        c.Label.String = 'Size of Largest Component';
    end
    %}
    
    degRatio_r = transpose(reshape(degRatio_r(:,1),[length(unique(pData_r(:,2))), length(unique(pData_r(:,1)))]));
    %{
    figure(1)
    subplot(2,3,r_ind)
    imagesc(degRatio_r)
    colormap('gray')
    c = colorbar
    caxis([0 1])
    title(['Conn. Radius =',num2str(uniqueRadi(r_ind))])
    set(gca,'XTick','')
    set(gca,'YTick','')
    
    if(r_ind == 4)
        set(gca,'XTick',1:2:9,'XTickLabel',unique_ihm(1:2:9))
        set(gca,'YTick',1:2:9,'YTickLabel',unique_ohm(1:2:9))
        xlabel('Incoming Half-Maximum')
        ylabel('Outgoing Half-Maximum')
        c.Label.String = 'Outgoing Degree Ratio';
    end
    %}
    
    count = 1;
    
%     for ii = 1:length(unique_ohm)
%        
%         find_ohm = pData_r(:,1) == unique_ohm(ii);
%         
%         for jj = 1:length(unique_ihm)
%             
%             find_ihm = pData_r(:,2) == unique_ihm(ii);
%             
%             findData = logical(find_ohm.*find_ihm);
%             
%             freqData_i = freqData_r{findData};
%             inFull_i = inFull_r{findData};
%             outFull_i = outFull_r{findData};
%             inDegDist_i = inDegDist_r{findData};
%             outDegDist_i = outDegDist_r{findData};
%             
%             freqData_outDeg_pairs = [outFull_i(freqData_i(:,1),2) freqData_i(:,3:4)];
%             freqData_inDeg_pairs = [inFull_i(freqData_i(:,1),2) freqData_i(:,3:4)];
%             
%             figure(r_ind)
%             subplot(10, 10, count)
%             %hold on
%             errorbar(freqData_inDeg_pairs(:,1), freqData_inDeg_pairs(:,2), freqData_inDeg_pairs(:,3),'.','Color',[0.85 0.33 0.1])
%             %errorbar(freqData_outDeg_pairs(:,1), freqData_outDeg_pairs(:,2), freqData_outDeg_pairs(:,3),'.','Color',[0 0.45 0.74])
%             
%             if(ii == 10 && jj == 1)
%                 xlabel('Degree')
%                 ylabel({['Out H.M = ',num2str(unique_ohm(ii))],'Frequency [Hz]'})
%             elseif(jj == 1)
%                 ylabel({['Out H.M = ',num2str(unique_ohm(ii))],''})
%             end
%             
%             if(ii == 1)
%                 title(['In H.M. = ',num2str(unique_ihm(jj))])
%             end
%             %}
%             
%             
%             %{
%             figure(r_ind)
%             subplot(10,10, count)
%             hold on
%             yyaxis left
%             plot(outDegDist_i(:,1), outDegDist_i(:,2),'.-')
% 
%             if(ii == 10 && jj == 1)
%                 xlabel('Degree')
%                 ylabel({['Out H.M = ',num2str(unique_ohm(ii))],'P(Degree)'})
%             elseif(jj == 1)
%                 ylabel({['Out H.M = ',num2str(unique_ohm(ii))],''})
%             end
%             
%             yyaxis right
%             plot(inDegDist_i(:,1), inDegDist_i(:,2),'.-')
%             hold off
%             
%             if(ii == 1)
%                 title(['In H.M. = ',num2str(unique_ihm(jj))])
%             end
%             %}   
%             
%             count = count + 1;
%             a = true;
%             
%         end
%         
%     end
    
    a = true;
    
end