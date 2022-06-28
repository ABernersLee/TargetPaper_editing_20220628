function plot_PV_Imagesc(corr_repeat_trial_session,all_labels,dirs,itype)

thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'PV_Imagesc/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

PlotPV = NaN(3,3,max(corr_repeat_trial_session(:,7)));

for irepeat = 1:3
%             irepeat2 = irepeat; %if irepeat>3; irepeat2 = 3; end
    
    if irepeat<3
        ind = corr_repeat_trial_session(:,5)==irepeat;
    else
        ind = corr_repeat_trial_session(:,5)>2;
    end
    
    oths = [1 find(~ismember(1:3,irepeat))+1];
    for isession = 1:max(corr_repeat_trial_session(:,7)) 
        mice = unique(corr_repeat_trial_session(corr_repeat_trial_session(:,7)==isession,8));
        
        for iothers = 1:length(oths)  
            
            PVmice = NaN(length(mice),1);
            for imouse = 1:length(mice)
                PVmice(imouse) = mean(corr_repeat_trial_session(ind & corr_repeat_trial_session(:,7)==isession & corr_repeat_trial_session(:,8)==mice(imouse),oths(iothers)),'omitnan');
            end
            if iothers==1
%                         PlotPV(irepeat,irepeat,isession) = nanmean(corr_repeat_trial_session(ind & corr_repeat_trial_session(:,7)==isession,oths(iothers)));                                                
                PlotPV(irepeat,irepeat,isession) = mean(PVmice);
            else
                PlotPV(irepeat,oths(iothers)-1,isession) = mean(PVmice);
            end
        end
    end                    
end

mn = mean(min(min(PlotPV)),'omitnan');
mx = mean(max(max(PlotPV)),'omitnan');

figure; hold on;
sz = floor(size(PlotPV,3)/3);
%     suptitle(['Early, middle and late sessions (n = ' num2str(size(PlotPV,3)) ', ' num2str(sz) ')'])
for isession = 1:3
    subplot(1,3,isession); hold on;

    set(gcf,'Position',[1706         649        1471         280])
    imagesc(mean(PlotPV(:,:,(1:sz)+(isession-1)*sz),3,'omitnan'))
    axis xy tight
    caxis([mn mx])       
    set(gca,'xtick',1:3,'xticklabel',all_labels.ittlabs(2:4))
    set(gca,'ytick',1:3,'yticklabel',all_labels.ittlabs(2:4))
    xlabel('Trial type')
    ylabel('Correlation with other trials')
    colorbar('Location','eastoutside')
    set(gca,'FontSize',8)
end        
helper_saveandclosefig([thisdir all_labels.typelab{itype} '_Imagesc_Over3SectionsOf' num2str(sz) all_labels.addon])

    