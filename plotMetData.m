clearvars
close all

load('metaData_2.mat');
load('metaData_2-degreesFull.mat');


unqOG_cutoff = unique(paramData(:,3));

for og_ind = 337:length(unqOG_cutoff)
   
    findOG = paramData(:,3) == unqOG_cutoff(og_ind);
    
    pData_og = paramData(findOG,[1 4]);
    odr_og = outDegreeRatio(findOG,:);
    odr_og_cv = odr_og(:,2)./odr_og(:,1);
    
    odr_og_r = flipud(transpose(reshape(odr_og(:,1),[length(unique(pData_og(:,2))), length(unique(pData_og(:,1)))])));
    %test = flipud(transpose(reshape(pData_og(:,1),[length(unique(pData_og(:,2))), length(unique(pData_og(:,1)))])));
    %{
    subplot(3,4,og_ind)
    imagesc(odr_og_r)
    colormap('gray')
    c = colorbar
    caxis([0 1])
    title(['Half-maximum =',num2str(unqOG_cutoff(og_ind))])
    set(gca,'XTick','')
    set(gca,'YTick','')
    
    if(og_ind == 9)
        set(gca,'XTick',1:2:9,'XTickLabel',unique_ig(1:2:9))
        set(gca,'YTick',1:length(unique(pData_og(:,1))),'YTickLabel',flipud(unique(pData_og(:,1))))
        xlabel('Incoming Half-Maximum')
        ylabel('Connection Radius')
        c.Label.String = 'Outgoing Degree Ratio';
    end
    %}
    
    unique_ig = unique(pData_og(:,2));
    
    lcs_og = cellfun('length', largestComponentNodes(findOG));
    lcs_og_r = flipud(transpose(reshape(lcs_og, [length(unique(pData_og(:,2))), length(unique(pData_og(:,1)))])));
    
    
    subplot(3,4,og_ind)
    imagesc(lcs_og_r./500)
    colormap('gray')
    c = colorbar;
    caxis([0 1])
    title(['Half-maximum =',num2str(unqOG_cutoff(og_ind))])
    set(gca,'XTick','')
    set(gca,'YTick','')
    
    if(og_ind == 9)
        set(gca,'XTick',1:2:9,'XTickLabel',unique_ig(1:2:9))
        set(gca,'YTick',1:length(unique(pData_og(:,1))),'YTickLabel',flipud(unique(pData_og(:,1))))
        xlabel('Incoming Half-Maximum')
        ylabel('Connection Radius')
        c.Label.String = 'Size of Largest Component';
    end
    %}
%     
%     freqData_og = freqData(findOG);
%     outDegDist_og = outDegreeDistribution(findOG);
%     inDegDist_og = inDegreeDistribution(findOG);
%     uniqueRange = unique(pData_og(:,1));
%     inFull_og = inDegree_full(findOG);
%     outFull_og = outDegree_full(findOG);
%     
%     maxDegree = 0;
%     
%     figure(og_ind)
%     
%     hAx=gobjects(length(uniqueRange)*length(unique_ig),1);
%     count = 1;
%     for r_ind = 1:length(uniqueRange)
%         
%         findRadi = pData_og(:,1)==uniqueRange(r_ind);
%         
%         id_radi = pData_og(findRadi,2);
%         
%         freqData_radi = freqData_og(findRadi);
%         outDegDist_radi = outDegDist_og(findRadi);
%         inDegDist_radi = inDegDist_og(findRadi);
%         inFull_radi = inFull_og(findRadi);
%         outFull_radi = outFull_og(findRadi);
%         
%         for id_ind = 1:length(id_radi)
%            
%             %{
%             freqData_i = freqData_radi{id_ind};
%             
%             freqData_outDeg_pairs = [outFull_radi{id_ind}(freqData_i(:,1),2) freqData_i(:,3:4)];
%             freqData_inDeg_pairs = [inFull_radi{id_ind}(freqData_i(:,1),2) freqData_i(:,3:4)];
%             
%             subplot(length(uniqueRange), length(id_radi), length(id_radi)*(r_ind-1)+id_ind)
%             hold on
%             errorbar(freqData_outDeg_pairs(:,1), freqData_outDeg_pairs(:,2),freqData_outDeg_pairs(:,3),'.','Color',[0.0 0.45 0.74])
%             errorbar(freqData_inDeg_pairs(:,1), freqData_inDeg_pairs(:,2),freqData_inDeg_pairs(:,3),'.','Color',[0.85 0.33 0.1])
%             hold off
%             
%             if(id_ind == 1)
%                 ylabel({['Conn. Radius = ',num2str(uniqueRange(r_ind))],''})
%                 if(r_ind==length(uniqueRange))
%                     ylabel({['Conn. Radius = ',num2str(uniqueRange(r_ind))],'Frquency [Hz]'})
%                     xlabel('Degree')
%                 end
%             end
%             
%             if(r_ind == 1)
%                 title(['In Half-Max. = ',num2str(round(id_radi(id_ind)/2)),' Hz'])
%             end
%             
%             a = true;
%             %}
%             
%             
%             %This commented block plots the degree distributions
%             hAx(count) = subplot(length(uniqueRange), length(id_radi), length(id_radi)*(r_ind-1)+id_ind);
%             hold on
%             yyaxis left
%             plot(outDegDist_radi{id_ind}(:,1),outDegDist_radi{id_ind}(:,2),'.-')
%             %ylim([0 1])
%             
%             if(id_ind == 1)
%                 ylabel({['Conn. Radius = ',num2str(uniqueRange(r_ind))],''})
%                 if(r_ind==length(uniqueRange))
%                     ylabel({['Conn. Radius = ',num2str(uniqueRange(r_ind))],'P(Degree)'})
%                     xlabel('Degree')
%                 end
%             end
%             
%             yyaxis right
%             plot(inDegDist_radi{id_ind}(:,1), inDegDist_radi{id_ind}(:,2),'.-')
%             %ylim([0 1])
%             hold off
%             %ylim([0 1])
%             
%             if(max([max(outDegDist_radi{id_ind}(:,1)) max(inDegDist_radi{id_ind}(:,1))]) > maxDegree)
%                 maxDegree = max([max(outDegDist_radi{id_ind}(:,1)) max(inDegDist_radi{id_ind}(:,1))]);
%             end
%             
%             if(r_ind == 1)
%                 title(['In Half-Max. = ',num2str(round(id_radi(id_ind)/2)),' Hz'])
%             end
%             %}
%             
%             count = count + 1;
%             
%         end
%         
%         a = true;
%         
%     end
%     
%     figure(og_ind)
%     %xlim(hAx,[-5 maxDegree+5])
%     
    a = true;
    
end