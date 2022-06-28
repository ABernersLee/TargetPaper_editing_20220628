function plot_Selectivity_AcrossOdors(selectivity_dat,selectivity_params,dirs) 

sp = selectivity_params;
sd = selectivity_dat;

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'AcrossOdors/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

%%%%% looking across the different probes and repeats to see if
%%%%% different numbers of neurons are significantly modulated to each
for iana = 3
    if iana==1
        P1 = sd.is_probe2; P2 = sd.session_probe2;
        R1 = sd.is_repeat2; R2 = sd.session_repeat2;
    elseif iana==3
        P1 = sd.is_probeMean2; P2 = sd.session_probeMean2;
        R1 = sd.is_repeatMean2; R2 = sd.session_repeatMean2;
          
        figure; hold on                        
        for ii = 1:2
            for imouse = 1:size(R1,1)   
                plot(ii+.2*imouse+rand(size(R1,2),1)/20,nansum(R1(imouse,:,:,ii),3)./nansum(R2(imouse,:,:,ii),3),'ok')                                
            end                
            errorbar(ii+.5,nanmean(nanmean(nansum(R1(:,:,:,ii),3)./nansum(R2(:,:,:,ii),3),2)),nanstd(nanmean(nansum(R1(:,:,:,ii),3)./nansum(R2(:,:,:,ii),3),2),1)./sqrt(size(R1,1)),'k','LineWidth',2)                
            yl = get(gca,'ylim');     
            text(ii+.5,yl(2)*.9,num2str(round(nanmean(nanmean(nansum(R1(:,:,:,ii),3)./nansum(R2(:,:,:,ii),3),2)),2,'significant')))                
        end                
        ylabel('Proportion of neurons significant')
        set(gca,'xtick',[1.5:2.5],'xticklabel',{'NT Repeat 1';'NT Repeat 2'})
        xlim([1 3])    
        
        helper_saveandclosefig([thisdir sp.ana_labs{iana} '_' sp.rep_probe{2}])
        
        figure; hold on                        
        for ii = 1:3
            for imouse = 1:size(P1,1)   
                plot(ii+.2*imouse+rand(size(P1,2),1)/20,nansum(P1(imouse,:,:,ii),3)./nansum(P2(imouse,:,:,ii),3),'ok')                                
            end                
            errorbar(ii+.5,nanmean(nanmean(nansum(P1(:,:,:,ii),3)./nansum(P2(:,:,:,ii),3),2)),nanstd(nanmean(nansum(P1(:,:,:,ii),3)./nansum(P2(:,:,:,ii),3),2),1)./sqrt(size(P1,1)),'k','LineWidth',2)                
            yl = get(gca,'ylim');     
            text(ii+.5,yl(2)*.9,num2str(round(nanmean(nanmean(nansum(P1(:,:,:,ii),3)./nansum(P2(:,:,:,ii),3),2)),2,'significant')))                
        end                
        ylabel('Proportion of neurons significant')
        set(gca,'xtick',[1.5:3.5],'xticklabel',{'Probe 1';'Probe 2';'Probe 3'})
        xlim([1 4])                    
        
        helper_saveandclosefig([thisdir sp.ana_labs{iana} '_' sp.rep_probe{1}])
        
    end
end

