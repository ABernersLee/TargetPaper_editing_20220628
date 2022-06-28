function plot_Selectivity_UpDownDuration(singledat,selectivity_dat,selectivity_params,dirs,figs_to_run)

od = singledat.toplotOdorLick(:,1:6000);

sp = selectivity_params;
sd = selectivity_dat;

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'UpDownDuration/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

     %%%%% Timing of when neurons come on and off
    titlab = {'Up, under ';'Up, >= ';'Down, all ';'Down, under ';'Down >= '};
    labslabs = {'odor-on up','odor-on down';'target','nontarget';'Up - target','Up - nontarget';'Up - target','Up - probe';'Up - probe','Up - nontarget';'Up - target','Up - repeats';'Up - repeat','Up - nontarget'};
    labind = [1 1 1;3 3 3;1 1 1;5:7;8:9 9];        
    for iana = [1 3]
        for itimes = 2
            if itimes==1
                times = sp.odorwindowlabel(:,1)>=0 & sp.odorwindowlabel(:,2)<=3000;
            else
                times = sp.odorwindowlabel(:,1)>=0 & sp.odorwindowlabel(:,2)<=2000;
            end
            timeslab = mean(sp.odorwindowlabel(times,1),2);

            for idat = 1:7
                 
                
                if sum(ismember(figs_to_run,[0 2 4]))==0 && idat==1
                    continue
                end

                if sum(ismember(figs_to_run,[0 4]))==0 && idat==3
                    continue
                end

                if sum(ismember(figs_to_run,0))==0 && (idat==2 || idat>3)
                    continue
                end

                odorbins = squeeze(singledat.odorsel(:,times,1,:));
            
                isUp = sd.all_Neurons(:,(idat-1)*2+3);
                isDown = sd.all_Neurons(:,(idat-1)*2+4);
                
                if iana==2
                    isUp = isUp & sd.neuron_info(:,7)==1;                    
                    isDown = isDown & sd.neuron_info(:,7)==1;
                elseif iana==3
                    if idat<3
                        isUp = sd.all_MeanSes(:,:,:,(idat-1)*2+3);
                        isDown = sd.all_MeanSes(:,:,:,(idat-1)*2+4);
                    elseif idat==3                        
                        isUp = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,5);
                        isDown = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,6);
                    elseif idat == 4                        
                        isUp = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,5);
                        isDown = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,10) & sd.all_MeanSes(:,:,:,9);
                    elseif idat == 5                        
                        isUp = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,10) & sd.all_MeanSes(:,:,:,9);                        
                        isDown = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,6);
                    elseif idat == 6
                        isUp = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,5);
                        isDown = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,12) & sd.all_MeanSes(:,:,:,11);                        
                    elseif idat == 7                        
                        isUp = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,12) & sd.all_MeanSes(:,:,:,11);
                        isDown = sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,6);
                    end
                    upCount = NaN(size(isUp,1),size(isUp,2),size(odorbins,2));
                    downCount = upCount;
                    upHist = NaN(size(isUp,1),size(isUp,2),size(odorbins,2));
                    downHist = upHist;
                    upCountOne = NaN(size(isUp));
                    for imouse = 1:4
                       [x,y] = find(squeeze(isUp(imouse,:,:)));
                       [x2,y2] = find(squeeze(isDown(imouse,:,:)));
                       xx = unique([x;x2]);
                       for isession = 1:length(xx)
                           indUp = sd.neuron_info(:,2)==imouse & sd.neuron_info(:,3)==xx(isession) & ismember(sd.neuron_info(:,6),y(x==xx(isession)));
                           indDown = sd.neuron_info(:,2)==imouse & sd.neuron_info(:,3)==xx(isession) & ismember(sd.neuron_info(:,6),y2(x2==xx(isession)));

                           if idat<3
                               upCountOne(imouse,xx(isession),y(xx(isession)==x)) = sum(odorbins(indUp,:,1)>0 & odorbins(indUp,:,2)<sp.p_cutoff,2);
                                   upCount(imouse,xx(isession),:) = histcounts(sum(odorbins(indUp,:,1)>0 & odorbins(indUp,:,2)<sp.p_cutoff,2),1:size(odorbins,2)+1)./sum(indUp);                           
                                   downCount(imouse,xx(isession),:) = histcounts(sum(odorbins(indDown,:,1)<0 & odorbins(indDown,:,2)<sp.p_cutoff,2),1:size(odorbins,2)+1)./sum(indDown);  
                                   upHist(imouse,xx(isession),:) = sum(odorbins(indUp,:,1)>0 & odorbins(indUp,:,2)<sp.p_cutoff,1)./sum(indUp);                           
                                   downHist(imouse,xx(isession),:) = sum(odorbins(indDown,:,1)<0 & odorbins(indDown,:,2)<sp.p_cutoff,1)./sum(indDown);
                           else                              
                               upCount(imouse,xx(isession),:) = histcounts(sum(odorbins(indUp,:,1)>0 & odorbins(indUp,:,2)<sp.p_cutoff,2),1:size(odorbins,2)+1)./sum(indUp);                           
                               downCount(imouse,xx(isession),:) = histcounts(sum(odorbins(indDown,:,1)>0 & odorbins(indDown,:,2)<sp.p_cutoff,2),1:size(odorbins,2)+1)./sum(indDown);  
                               upHist(imouse,xx(isession),:) = sum(odorbins(indUp,:,1)>0 & odorbins(indUp,:,2)<sp.p_cutoff,1)./sum(indUp);                           
                               downHist(imouse,xx(isession),:) = sum(odorbins(indDown,:,1)>0 & odorbins(indDown,:,2)<sp.p_cutoff,1)./sum(indDown);
                           end
                       end
                    end                        
                end
                
                        
                if iana<3

                    upCount = sum(odorbins(isUp,:,1)>0 & odorbins(isUp,:,2)<(sp.p_cutoff/size(odorbins,2)),2);
                    downCount = sum(odorbins(isDown,:,1)<0 & odorbins(isDown,:,2)<(sp.p_cutoff/size(odorbins,2)),2);
                    upHist = sum(odorbins(isUp,:,1)>0 & odorbins(isUp,:,2)<(sp.p_cutoff/size(odorbins,2)),1);
                    downHist = sum(odorbins(isDown,:,1)<0 & odorbins(isDown,:,2)<(sp.p_cutoff/size(odorbins,2)),1);
                
                    if sum(ismember(figs_to_run,0))>0
                
                        figure; hold on
                        subplot(2,2,1); hold on
                        histogram(upCount(upCount>0),1:size(odorbins,2)-1)
                        set(gca,'xtick',1:size(odorbins,2)-1)
                        ylabel(['Number of ' labslabs{idat,1} ' cells'])
                        xlabel('Duration of modulation, number sig. up bins (200ms each)')
    
                        subplot(2,2,3); hold on
                        histogram(downCount(downCount>0),1:size(odorbins,2)-1)
                        set(gca,'xtick',1:2:size(odorbins,2)-1)
                        ylabel(['Number of ' labslabs{idat,2} ' cells'])
                        xlabel('Duration of modulation, number sig. down bins (200ms each)')
    
    
                        subplot(2,2,2); hold on
                        bar(upHist)
                        set(gca,'xtick',1:2:length(upHist),'xticklabel',timeslab(1:2:end))        
                        ylabel(['Number of ' labslabs{idat,1} ' cells sig per bin'])
                        xlabel('Bin start (ms)')
    
                        subplot(2,2,4); hold on
                        bar(downHist)
                        set(gca,'xtick',1:2:length(upHist),'xticklabel',timeslab(1:2:end))
                        ylabel(['Number of ' labslabs{idat,2} ' cells sig per bin'])
                        xlabel('Bin start (ms)')
    
                        suptitle([num2str(sum(times)*.2) ' Seconds, ' sp.ana_labs{iana} ', type ' num2str(idat)])
    %                     set(gcf,'Position',[2036         220         853         849])
                    end
                else
                    if (sum(ismember(figs_to_run,[0 2]))>0) || (sum(ismember(figs_to_run,4))>0 && idat==3)

                    figure; hold on
                    subplot(2,2,1); hold on           
                    for imouse = 1:4
%                         plot(repmat([1:size(odorbins,2)]',[1 size(upCount,2)])+rand(size(odorbins,2),size(upCount,2))/5,squeeze(upCount(imouse,:,:))','.','color',[.5 .5 .5])
                        plot(1:size(odorbins,2)-1,nanmean(squeeze(upCount(imouse,:,1:end-1)),1),'-o','color',[.5 .5 .5])        
                        errorbar([1:size(odorbins,2)-1],nanmean(squeeze(upCount(imouse,:,1:end-1)),1)',nanstd(squeeze(upCount(imouse,:,1:end-1)),1)'./squeeze(sqrt(sum(~isnan(upCount(imouse,:,1:end-1)),2))),'color',[.5 .5 .5])     
                    end
                    bar(1:size(odorbins,2)-1,squeeze(nanmean(nanmean(upCount(:,:,1:end-1),1),2)),'FaceColor','none','EdgeColor','k')                    
                    errorbar(1:size(odorbins,2)-1,squeeze(nanmean(nanmean(upCount(:,:,1:end-1),1),2)),nanstd(squeeze(nanmean(upCount(:,:,1:end-1),2)),1)./sqrt(size(upCount(:,:,1:end-1),1)),'k','Linewidth',2)                    
                    set(gca,'xtick',1:2:size(odorbins,2)+1)
                    ylabel(['Number of ' labslabs{idat,1} ' cells'])
                    xlabel('Duration of modulation, number sig. up bins (200ms each)')

                    subplot(2,2,3); hold on     
                    for imouse = 1:4
                        plot(1:size(odorbins,2)-1,nanmean(squeeze(downCount(imouse,:,1:end-1)),1),'-o','color',[.5 .5 .5])        
                        errorbar([1:size(odorbins,2)-1],nanmean(squeeze(downCount(imouse,:,1:end-1)),1)',nanstd(squeeze(downCount(imouse,:,1:end-1)),1)'./squeeze(sqrt(sum(~isnan(downCount(imouse,:,1:end-1)),2))),'color',[.5 .5 .5])     
                    end
                    bar(1:size(odorbins,2)-1,squeeze(nanmean(nanmean(downCount(:,:,1:end-1),1),2)),'FaceColor','none','EdgeColor','k')                    
                    errorbar(1:size(odorbins,2)-1,squeeze(nanmean(nanmean(downCount(:,:,1:end-1),1),2)),nanstd(squeeze(nanmean(downCount(:,:,1:end-1),2)),1)./sqrt(size(downCount(:,:,1:end-1),1)),'k','Linewidth',2)                   
                    set(gca,'xtick',1:2:size(odorbins,2)+1)
                    ylabel(['Number of ' labslabs{idat,2} ' cells'])
                    xlabel('Duration of modulation, number sig. down bins (200ms each)')
                    
                    dat1 = upCount(:,:,1:end-1);
                    dat2 = downCount(:,:,1:end-1);
                    an2mouse = repmat([1:size(dat1,1)]',[1 size(dat1,2) size(dat1,3)]);
                    an2session = permute(repmat([1:size(dat1,2)]',[1 size(dat1,1) size(dat1,3)]),[2 1 3]);                    
                    an2timebin = permute(repmat([1:size(dat1,3)]',[1 size(dat1,1) size(dat1,2)]),[2 3 1]);
                    [pLeft,qLeft,~] = anovan([dat1(:);dat2(:)],[[ones(size(dat1(:)));2*ones(size(dat2(:)))]...
                        [an2mouse(:);an2mouse(:)] [an2session(:);an2session(:)]...
                        [an2timebin(:);an2timebin(:)]],'varnames',{'Group','Mouse','Session','TimeBin'},...
                        'continuous',[3 4],'model',[eye(4,4);[1 0 0 1]],'display','off');
%                         'continuous',[3 4],'model','interaction','display','off');
                    title(['Anova groupXtime inter. p = ' num2str(round(pLeft(5),2,'significant'))])
                    
                    subplot(2,2,2); hold on
                    for imouse = 1:4
                        plot(1:size(odorbins,2),nanmean(squeeze(upHist(imouse,:,:)),1),'-o','color',[.5 .5 .5])        
                        errorbar([1:size(odorbins,2)],nanmean(squeeze(upHist(imouse,:,:)),1)',nanstd(squeeze(upHist(imouse,:,:)),1)'./squeeze(sqrt(sum(~isnan(upHist(imouse,:,:)),2))),'color',[.5 .5 .5])     
                    end                    
                    bar(1:size(odorbins,2),squeeze(nanmean(nanmean(upHist,1),2)),'FaceColor','none','EdgeColor','k')                    
                    errorbar(1:size(odorbins,2),squeeze(nanmean(nanmean(upHist,1),2)),nanstd(squeeze(nanmean(upHist,2)),1)./sqrt(size(upHist,1)),'k','Linewidth',2)                        
                    set(gca,'xtick',1:length(upHist),'xticklabel',timeslab(1:end))        
                    ylabel(['Number of ' labslabs{idat,1} ' cells sig per bin'])
                    xlabel('Bin start (ms)')

                    subplot(2,2,4); hold on
                    for imouse = 1:4
                        plot(1:size(odorbins,2),nanmean(squeeze(downHist(imouse,:,:)),1),'-o','color',[.5 .5 .5])        
                        errorbar([1:size(odorbins,2)],nanmean(squeeze(downHist(imouse,:,:)),1)',nanstd(squeeze(downHist(imouse,:,:)),1)'./squeeze(sqrt(sum(~isnan(downHist(imouse,:,:)),2))),'color',[.5 .5 .5])     
                    end                    
                    bar(1:size(odorbins,2),squeeze(nanmean(nanmean(downHist,1),2)),'FaceColor','none','EdgeColor','k')                    
                    errorbar(1:size(odorbins,2),squeeze(nanmean(nanmean(downHist,1),2)),nanstd(squeeze(nanmean(downHist,2)),1)./sqrt(size(downHist,1)),'k','Linewidth',2)    
                    set(gca,'xtick',1:length(downHist),'xticklabel',timeslab(1:end))
                    ylabel(['Number of ' labslabs{idat,2} ' cells sig per bin'])
                    xlabel('Bin start (ms)')
                    
                    dat1 = upHist;
                    dat2 = downHist;
                    an2mouse = repmat([1:size(dat1,1)]',[1 size(dat1,2) size(dat1,3)]);
                    an2session = permute(repmat([1:size(dat1,2)]',[1 size(dat1,1) size(dat1,3)]),[2 1 3]);                    
                    an2timebin = permute(repmat([1:size(dat1,3)]',[1 size(dat1,1) size(dat1,2)]),[2 3 1]);
                    [pRight,qRight,~] = anovan([dat1(:);dat2(:)],[[ones(size(dat1(:)));2*ones(size(dat2(:)))]...
                        [an2mouse(:);an2mouse(:)] [an2session(:);an2session(:)]...
                        [an2timebin(:);an2timebin(:)]],'varnames',{'Group','Mouse','Session','TimeBin'},...
                        'continuous',[3 4],'model',[eye(4,4);[1 0 0 1]],'display','off');
%                     'continuous',[3 4],'model','interaction','display','off');
                    title(['Anova groupXtime inter. p = ' num2str(round(pRight(5),2,'significant'))])

                    suptitle([num2str(sum(times)*.2) ' Seconds, ' sp.ana_labs{iana} ', type ' num2str(idat)])
%                     set(gcf,'Position',[2036         220         853         849])
                    end
                end
                
                if  (iana==3 && sum(ismember(figs_to_run,[0 2]))>0) || (sum(ismember(figs_to_run,4))>0 && idat==3 && iana==3) || (sum(ismember(figs_to_run,0))>0)                   
                    set(gcf,'Position',[680   324   958   654])
                    helper_saveandclosefig([thisdir 'DistDurationSig_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_type' num2str(idat) '_' sp.ana_labs{iana}])
                    
                    if iana>=3
                                            
%                         writecell(qRight,[thisdir 'DistDurationSig_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_type' num2str(idat) '_' sp.ana_labs{iana} '_ANOVAN_right.txt'],'Delimiter','tab')                        
                        writecell(qRight,[thisdir 'DistDurationSig_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_type' num2str(idat) '_' sp.ana_labs{iana} '_ANOVAN_right.xlsx'])   
                        writecell(qLeft,[thisdir 'DistDurationSig_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_type' num2str(idat) '_' sp.ana_labs{iana} '_ANOVAN_left.xlsx'])   
                        
                    end
                else
                    close all
                end
                

                
                if itimes==2 && idat==1 && iana==3
                    

                   if sum(ismember(figs_to_run,[0 4]))>0 
                        %proportions target or nontarget
                        isUpTargetNT = sd.all_MeanSes(:,:,:,3) & (sd.all_MeanSes(:,:,:,5) | sd.all_MeanSes(:,:,:,6)) & ~sp.dem_num;             
                        isUpNot = sd.all_MeanSes(:,:,:,3) & ~sd.all_MeanSes(:,:,:,5) & ~sd.all_MeanSes(:,:,:,6) & ~sp.dem_num;
                        isDownTargetNT = sd.all_MeanSes(:,:,:,4) & (sd.all_MeanSes(:,:,:,5) | sd.all_MeanSes(:,:,:,6)) & ~sp.dem_num;           
                        isDownNot = sd.all_MeanSes(:,:,:,4) & ~sd.all_MeanSes(:,:,:,5) & ~sd.all_MeanSes(:,:,:,6) & ~sp.dem_num;
                        
                        isBriefTarget = upCountOne<sp.medcut & sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,5) & ~sd.all_MeanSes(:,:,:,6) & ~sp.dem_num;              
                        isBriefNT = upCountOne<sp.medcut & sd.all_MeanSes(:,:,:,3) & ~sd.all_MeanSes(:,:,:,5) & sd.all_MeanSes(:,:,:,6) & ~sp.dem_num;
                        isLongTarget = upCountOne>=sp.medcut & sd.all_MeanSes(:,:,:,3) & sd.all_MeanSes(:,:,:,5) & ~sd.all_MeanSes(:,:,:,6) & ~sp.dem_num;              
                        isLongNT = upCountOne>=sp.medcut & sd.all_MeanSes(:,:,:,3) & ~sd.all_MeanSes(:,:,:,5) & sd.all_MeanSes(:,:,:,6) & ~sp.dem_num;
                        
                        figure; hold on
                        dat1 = nansum(isUpTargetNT,3)./nansum((isUpNot | isUpTargetNT),3);
                        dat2 = nansum(isDownTargetNT,3)./nansum((isDownNot | isDownTargetNT),3);
                        plot(1:4,dat1,'ko')                    
                        plot(6:9,dat2,'ro')
                        dat = nansum(nansum(isUpTargetNT,3),2)./nansum(nansum(isUpNot | isUpTargetNT,3),2);
                        errorbar(2.5,nanmean(dat),nanstd(dat)./sqrt(length(dat)),'k','LineWidth',2)                    
                        dat = nansum(nansum(isDownTargetNT,3),2)./nansum(nansum(isDownNot | isDownTargetNT,3),2);
                        errorbar(7.5,nanmean(dat),nanstd(dat)./sqrt(length(dat)),'r','LineWidth',2)
                        set(gca,'xtick',[2.5 7.5],'xticklabel',{'Up';'Down'})
                        ylabel('Proportion of sig target/NT mod neurons')
                        xlim([0 10])
                        ylim([-.1 1.1]) 
                        
                        an2mouse = repmat([1:size(dat1,1)]',[1 size(dat1,2)]);
                        an2session = permute(repmat([1:size(dat1,2)]',[1 size(dat1,1)]),[2 1]);      
                        [p,tbl,~] = anovan([dat1(:);dat2(:)],[[ones(size(dat1(:)));2*ones(size(dat2(:)))]...
                            [an2mouse(:);an2mouse(:)] [an2session(:);an2session(:)]]...
                            ,'display','off','varnames',{'UpDown';'Mouse';'Session'},'continuous',3);
                        title(['Anova group p = ' num2str(round(p(1),2,'significant'))])
                        
                        
                       helper_saveandclosefig([thisdir 'PropSigUpDownTNT_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_type' num2str(idat) '_' sp.ana_labs{iana}])
                       writecell(tbl,[thisdir 'PropSigUpDownTNT_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_type' num2str(idat) '_' sp.ana_labs{iana} '_ANOVA.txt'],'Delimiter','tab')               
                       writecell(tbl,[thisdir 'PropSigUpDownTNT_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_type' num2str(idat) '_' sp.ana_labs{iana} '_ANOVA.xlsx'])               
                    end

                   if sum(ismember(figs_to_run,0))>0                             
                        figure; hold on
                        dat1 = nansum(isBriefTarget,3)./nansum((isBriefNT | isBriefTarget),3);
                        dat2 = nansum(isLongTarget,3)./nansum((isLongNT | isLongTarget),3);
    %                     plot(1:4,dat1,'ko')                    
    %                     plot(6:9,dat2,'ro')
                        errorbar(1:4,nanmean(dat1,2),nanstd(dat1,[],2)./sqrt(sum(~isnan(dat1),2)),'k','LineStyle','none')                    
                        errorbar(6:9,nanmean(dat2,2),nanstd(dat2,[],2)./sqrt(sum(~isnan(dat2),2)),'r','LineStyle','none')
                        dat = nansum(nansum(isBriefTarget,3),2)./nansum(nansum(isBriefNT | isBriefTarget,3),2);
                        errorbar(2.5,nanmean(dat),nanstd(dat)./sqrt(length(dat)),'k','LineWidth',2)                    
                        dat = nansum(nansum(isLongTarget,3),2)./nansum(nansum(isLongNT | isLongTarget,3),2);
                        errorbar(7.5,nanmean(dat),nanstd(dat)./sqrt(length(dat)),'r','LineWidth',2)
                        set(gca,'xtick',[2.5 7.5],'xticklabel',{'Up Brief';'Up Long'})
                        ylabel('Proportion of Target preferring neurons of all sig mod')                    
                        xlim([0 10])
                        ylim([-.1 1.1]) 
                                          
                        an2mouse = repmat([1:size(dat1,1)]',[1 size(dat1,2)]);
                        an2session = permute(repmat([1:size(dat1,2)]',[1 size(dat1,1)]),[2 1]);      
                        [p,~,~] = anovan([dat1(:);dat2(:)],[[ones(size(dat1(:)));2*ones(size(dat2(:)))]...
                            [an2mouse(:);an2mouse(:)] [an2session(:);an2session(:)]]...
                            ,'display','off');
                        title(['Anova group p = ' num2str(round(p(1),2,'significant'))])
                        
                        helper_saveandclosefig([thisdir 'PropSigBriefLongTNT_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_type' num2str(idat) '_' sp.ana_labs{iana} '_medcut' num2str(sp.medcut)])
                   end
                end
                    
                
                if itimes==2 && idat==1 && iana==1 && sum(ismember(figs_to_run,[0 2]))>0 
                    fUp = find(isUp); fDown = find(isDown);                    
                    alldat2 = {fUp(upCount<sp.medcut & upCount>0),fUp(upCount>=sp.medcut),fDown,fDown(downCount<sp.medcut),fDown(downCount>=sp.medcut)};
                    figure; hold on                    
                    set(gcf,'Position',[ 1951         331         987        1033])
                    clear subp; yls = NaN(5,2);
                    for idat2 = 1:5                        
                        subplot(5,2,(idat2-1)*2+1); hold on
                        n = alldat2{idat2};                         
                        title([titlab{idat2} num2str(sp.medcut) ' bins, n=' num2str(length(n))])
                        if length(n)>1
                            dat = od(n,:);                           
                            [~,mm] = max(dat(:,3000:-1:1000),[],2);
                            [~,ord] = sort(mm);
                            
                            dat2 = zscore(dat,[],2);
                            datBS = dat2-mean(dat2(:,1:1000),2);
                            

%                             imagesc(dat(ord,:)./max(dat(ord,:),[],2))
                            % set(gca,'ytick',1:size(dat,1),'yticklabel',n(ord))
                            imagesc(dat(ord,:)./max(dat(ord,:),[],2))
                            set(gca,'ytick',[])
                            ylabel('Each cell')
                            xlabel('Seconds from odor-on')
                            set(gca,'xtick',1000:1000:6000,'xticklabel',0:5) 
                            axis tight ij
                            colorbar

                            subp(idat2) = subplot(5,2,(idat2-1)*2+2); hold on
                            tbs = 1:size(dat,2);
%                             m1 = mean(zscore(dat,[],2)); sem1 = std(zscore(dat,[],2))./sqrt(size(dat,1));
                            m1 = mean(datBS); sem1 = std(datBS)./sqrt(size(datBS,1));
                            m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
                            rev = m+sem; fwd = m-sem;
                            p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
                            p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
                            plot(mean(datBS,'omitnan'),'k');
                            set(gca,'xtick',1000:1000:6000,'xticklabel',0:5)
                            axis tight
                            yl = get(gca,'ylim');
%                             plot([1000 1000],yl,'r')
%                             plot([3000 3000],yl,'r')
                            ylabel('Avg z-scored FR')
                            xlabel('Seconds from odor-on')
                            yls(idat2,:) = yl;
                        end
                    end

                    set(gcf,'Position',[ 1951         331         987        1033])
                    for idat2 = 1:5
                        if idat2<3
                            subp(idat2).YLim = [min(yls(1:2,1)) max(yls(1:2,2))];
                        end
                        axes(subp(idat2))
                        plot([1000 1000],subp(idat2).YLim,'r')
                        plot([3000 3000],subp(idat2).YLim,'r')
                    end
                    helper_saveandclosefig([thisdir 'BriefLong_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_medcut' num2str(sp.medcut) '_type' num2str(idat) '_' sp.ana_labs{iana}])

                    figure; hold on

                    data1 = mean(od(isUp,1:500),2)*1000;
                    errorbar(1,mean(data1),std(data1)./sqrt(sum(~isnan(data1))),'k')
                    data2 = mean(od(isDown,1:500),2)*1000;
                    errorbar(2,mean(data2),std(data2)./sqrt(sum(~isnan(data2))),'k')
                    set(gca,'xlim',[.5 2.5])
                    p = ranksum(data1,data2);
                    ylabel('FR'); set(gca,'xtick',1:2,'xticklabel',...
                        {['Excitied (N = ' num2str(sum(isUp)) ')'],...
                        ['Inhibited (N = ' num2str(sum(isDown)) ')']})
                    title(['The 1 second before odor, RanksumP = ' num2str(p)])
                      helper_saveandclosefig([thisdir 'BriefLong_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_medcut' num2str(sp.medcut) '_type' num2str(idat) '_' sp.ana_labs{iana} '_other'])

                end
                
            end
        end
    end

  