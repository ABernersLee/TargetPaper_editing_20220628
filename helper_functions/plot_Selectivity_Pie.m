function plot_Selectivity_Pie(selectivity_dat,selectivity_params,dirs,figs_to_run)

sp = selectivity_params;
sd = selectivity_dat;

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'Pie/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end
plotlab2 = {'';'ExclLickNeurons'};

%%%%% Pie charts and errorbars for how many neurons are 
    % significant in certain things
    
for iana = [1 3] 
    if iana == 1
        islick = sd.all_Neurons(:,2);    
        isUp = sd.all_Neurons(:,3);            
        isDown = sd.all_Neurons(:,4);
        isOdorOff = sd.all_Neurons(:,13);
        PS = sd.all_Neurons(:,9);
        PA = sd.all_Neurons(:,10);
        RS = sd.all_Neurons(:,11);
        RA = sd.all_Neurons(:,12);
        TNT = sd.all_Neurons(:,5) | sd.all_Neurons(:,6);
        islickU = sd.islickup; islickD = sd.islickdown;
    elseif iana == 2            
        islick = reshape(sd.all_BestSes(:,:,2),[size(sd.all_BestSes,1)*size(sd.all_BestSes,2) 1]);
        isUp = reshape(sd.all_BestSes(:,:,3),[size(sd.all_BestSes,1)*size(sd.all_BestSes,2) 1]);
        isDown = reshape(sd.all_BestSes(:,:,4),[size(sd.all_BestSes,1)*size(sd.all_BestSes,2) 1]);          
        isOdorOff = reshape(sd.all_BestSes(:,:,13),[size(sd.all_BestSes,1)*size(sd.all_BestSes,2) 1]);            
        PS = reshape(sd.all_BestSes(:,:,9),[size(sd.all_BestSes,1)*size(sd.all_BestSes,2) 1]);
        PA = reshape(sd.all_BestSes(:,:,10),[size(sd.all_BestSes,1)*size(sd.all_BestSes,2) 1]);
        RS = reshape(sd.all_BestSes(:,:,11),[size(sd.all_BestSes,1)*size(sd.all_BestSes,2) 1]);
        RA = reshape(sd.all_BestSes(:,:,12),[size(sd.all_BestSes,1)*size(sd.all_BestSes,2) 1]);         
        TNT =reshape(sd.all_BestSes(:,:,5) | sd.all_BestSes(:,:,6) ,[size(sd.all_BestSes,1)*size(sd.all_BestSes,2) 1]);
    elseif iana==3            
        islick = sd.all_MeanSes(:,:,:,2);   
        isUp = sd.all_MeanSes(:,:,:,3);            
        isDown = sd.all_MeanSes(:,:,:,4);
        isOdorOff = sd.all_MeanSes(:,:,:,13);
        PS = sd.all_MeanSes(:,:,:,9);
        PA = sd.all_MeanSes(:,:,:,10);
        RS = sd.all_MeanSes(:,:,:,11);
        RA = sd.all_MeanSes(:,:,:,12);
        TNT = sd.all_MeanSes(:,:,:,5) | sd.all_MeanSes(:,:,:,6);
        islickU = sd.all_MeanSes(:,:,:,17);  islickD = sd.all_MeanSes(:,:,:,18);
    end

    if iana<3

        if sum(ismember(figs_to_run,0))>0
            figure; hold on
            subplot(3,2,1)
            pie([sum(~TNT), sum(TNT),],{'Not Modulated';'Target Modulated'})            
            subplot(3,2,2)
            pie([sum(~TNT), sum(TNT),])            
            subplot(3,2,3)
            pie([sum(PA & ~PS), sum(PA & PS),],{'Not Modulated';'Probe Modulated'})            
            subplot(3,2,4)
            pie([sum(PA & ~PS), sum(PA & PS),])            
            subplot(3,2,5)
            pie([sum(RA & ~RS), sum(RA & RS),],{'Not Modulated';'Repeat Modulated'})            
            subplot(3,2,6)
            pie([sum(RA & ~RS), sum(RA & RS),])
            helper_saveandclosefig([thisdir 'PieModulations_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana}])
        end
         
        if sum(ismember(figs_to_run,[0 2 15 14]))>0
            figure; hold on
            subplot(2,3,1)
            pie([sum(~isUp & ~isDown), sum(isUp), sum(isDown)],{'Not Modulated';'Odor-on up';'Odor-on down'})
            subplot(2,3,2)
            pie([sum(~isUp & ~isDown & ~islick),sum(islick), sum(isUp), sum(isDown)],{'Not Modulated';'Lick-on';'Odor-On up';'Odor-On down'})
            subplot(2,3,3)
            pie([sum(~isOdorOff & ~isDown & ~isUp), sum(isOdorOff & ~isDown & ~isUp),sum(isOdorOff & (isUp | isDown)), sum(~isOdorOff & isUp & ~isDown),sum(~isOdorOff & ~isUp & isDown)],...
                {'Not Modulated';'Odor-off only';'Odor-off and odor-on';'Odor-on up only';'Odor-on down only'})
            subplot(2,3,4)
            pie([sum(~isOdorOff & ~isDown & ~isUp), sum(isOdorOff & ~isDown & ~isUp),sum(isOdorOff & isUp & ~isDown),sum(isOdorOff & ~isUp & isDown), sum(~isOdorOff & isUp & ~isDown),sum(~isOdorOff & ~isUp & isDown)],...
                {'Not Modulated';'Odor-off only';'Odor-off and odor-on up';'Odor-off and odor-on down';'Odor-on up only';'Odor-on down only'})            
            subplot(2,3,5)
            pie([sum(~isUp & ~isDown & ~islick),sum(islickU), sum(islickD), sum(isUp), sum(isDown)],{'Not Modulated';'Lick-on up';'Lick-on down';'Odor-On up';'Odor-On down'})
            set(gcf,'Position',[1947         253        1017         711])
            helper_saveandclosefig([thisdir 'PieWithLabels_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana}])

            figure; hold on
            subplot(2,3,1)
            pie([sum(~isUp & ~isDown), sum(isUp), sum(isDown)],'%.3f%%')
            subplot(2,3,2)
            pie([sum(~isUp & ~isDown & ~islick),sum(islick), sum(isUp), sum(isDown)],'%.3f%%')
            subplot(2,3,3)
            pie([sum(~isOdorOff & ~isDown & ~isUp), sum(isOdorOff & ~isDown & ~isUp),sum(isOdorOff & (isUp | isDown)), sum(~isOdorOff & isUp & ~isDown),sum(~isOdorOff & ~isUp & isDown)],'%.3f%%')
            subplot(2,3,4)
            pie([sum(~isOdorOff & ~isDown & ~isUp), sum(isOdorOff & ~isDown & ~isUp),sum(isOdorOff & isUp & ~isDown),sum(isOdorOff & ~isUp & isDown), sum(~isOdorOff & isUp & ~isDown),sum(~isOdorOff & ~isUp & isDown)],'%.3f%%')
              subplot(2,3,5)
            pie([sum(~isUp & ~isDown & ~islick),sum(islickU), sum(islickD), sum(isUp), sum(isDown)],'%.3f%%')
            set(gcf,'Position',[1947         253        1017         711])
            helper_saveandclosefig([thisdir 'PieWithNumbers_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana}])
        end
         
    else
        

        if sum(ismember(figs_to_run,[0 4]))>0
            for iplot2 = 1:2        
                figure; hold on
                 dat = zeros(size(isUp,1),size(isUp,2),size(isUp,3),5);
                 if iplot2 == 1
                    IND = true(size(islick));
                 else
                    IND = ~islick;
                 end
                dat(:,:,:,1) = TNT & IND; 
                dat(:,:,:,2) = PS & IND;                 
                dat(:,:,:,3) = RS & IND; 
                dat(:,:,:,4) = PA & IND;                 
                dat(:,:,:,5) = RA & IND;     
                dat(repmat(sp.dem_num,[1 1 1 5])) = NaN;
                for imouse = 1:size(isUp,1)                                
                    plot(1+.2*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,1),3),'ok')
                    plot(2+.2*imouse+rand(size(isUp,2),1)/20,nansum(dat(imouse,:,:,2),3)./nansum(dat(imouse,:,:,4),3),'ok')                
                    plot(3+.2*imouse+rand(size(isUp,2),1)/20,nansum(dat(imouse,:,:,3),3)./nansum(dat(imouse,:,:,5),3),'ok')                                 
                end                                     
                for iplot = 1:3
                    if iplot==1
                        errorbar(iplot+.5,nanmean(nanmean(nanmean(dat(:,:,:,iplot),3),2)),std(nanmean(nanmean(dat(:,:,:,iplot),3),2),1)./sqrt(size(dat,1)),'k','LineWidth',2)
                    elseif iplot>1                    
                        errorbar(iplot+.5,nanmean(nanmean(nansum(dat(:,:,:,iplot),3)./nansum(dat(:,:,:,iplot+2),3),2)),nanstd(nanmean(nansum(dat(:,:,:,iplot),3)./nansum(dat(:,:,:,iplot+2),3),2),1)./sqrt(size(dat,1)),'k','LineWidth',2)
                    end
                end          
                ylabel('Proportion of neurons significant')
                 set(gca,'xtick',[1.5:3.5],'xticklabel',{'Target';'Probe';'NT Repeat'})
                 yl = get(gca,'ylim');
                for iplot = 1:3
                    if iplot==1
                        text(iplot+.5,yl(2)*.9,num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,iplot),3),2)),2,'significant')))
                    else
                        text(iplot+.5,yl(2)*.9,num2str(round(nanmean(nanmean(nansum(dat(:,:,:,iplot),3)./nansum(dat(:,:,:,iplot+2),3),2)),2,'significant')))
                    end
    
                end
                dat1 = nanmean(dat(:,:,:,1),3);
                dat2 = nansum(dat(:,:,:,2),3)./nansum(dat(:,:,:,4),3);
                dat3 = nansum(dat(:,:,:,3),3)./nansum(dat(:,:,:,5),3);
                an2mouse = repmat([1:size(dat1,1)]',[1 size(dat1,2)]);            
                an2session = repmat([1:size(dat1,2)]',[1 size(dat1,1)])';
                [p,tbl,stats] = anovan([dat1(:);dat2(:);dat3(:)],...
                    [[ones(size(dat1(:)));2*ones(size(dat2(:)));3*ones(size(dat3(:)))] [an2mouse(:);an2mouse(:);an2mouse(:)] [an2session(:);an2session(:);an2session(:)]],...
                    'display','off','varnames',{'TrialType';'Mouse';'Session'},'continuous',3);
                stats2 = multcompare(stats,'display','off','ctype','bonferroni');
                ps = stats2(:,end);
                title(['Anova with session, mouse and type. Type p = ' num2str(round(p(1),2,'significant'))])
                if ps(1)<.05
                    text(1.1,yl(2)*1.1,['Target v Probe: ' num2str(round(ps(1),2,'significant'))])
                end
                if ps(2)<.05
                    text(2.5,yl(2)*1.1,['Target v Repeats: ' num2str(round(ps(2),2,'significant'))])
                end            
                if ps(3)<.05
                    text(2.5,yl(2)*1,['Probe v Repeats: ' num2str(round(ps(3),2,'significant'))])
                end
                ylim([yl(1) yl(2)*1.2])
                helper_saveandclosefig([thisdir 'PieModulations_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} plotlab2{iplot2}])
                writecell(tbl,[thisdir 'PieModulations_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} plotlab2{iplot2} '_ANOVA.txt'],'Delimiter','tab') 
                writecell(num2cell(stats2),[thisdir 'PieModulations_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} plotlab2{iplot2} '_PostHoc.txt'],'Delimiter','tab')       
                writecell(tbl,[thisdir 'PieModulations_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} plotlab2{iplot2} '_ANOVA.xlsx']) 
                writecell(num2cell(stats2),[thisdir 'PieModulations_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} plotlab2{iplot2} '_PostHoc.xlsx'])       
            end
        end

        if sum(ismember(figs_to_run,[0 2]))>0
            figure; hold on;   
            dat = zeros(size(isUp,1),size(isUp,2),size(isUp,3),2);
            dat(:,:,:,1) = isUp; 
            dat(:,:,:,2) = isDown;     
            dat(repmat(sp.dem_num,[1 1 1 2])) = NaN;
            for imouse = 1:size(isUp,1)                
                plot(2+.1*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,1),3),'ok') 
                plot(3+.1*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,2),3),'ok')                 
            end            
            errorbar(2.5,nanmean(nanmean(nanmean(dat(:,:,:,1),3),2)),std(nanmean(nanmean(dat(:,:,:,1),3),2),1)./sqrt(size(dat(:,:,:,1),1)),'k','LineWidth',2)
            errorbar(3.5,nanmean(nanmean(nanmean(dat(:,:,:,2),3),2)),std(nanmean(nanmean(dat(:,:,:,2),3),2),1)./sqrt(size(dat(:,:,:,2),1)),'k','LineWidth',2)                                
            set(gca,'xtick',[2.25:3.25],'xticklabel',{'Odor-on up';'Odor-on down'})
            yl = get(gca,'ylim');
            text(2.5,yl(2)*.99,[num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,1),3),2)),2,'significant'))])           
            text(3.5,yl(2)*.99,[num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,2),3),2)),2,'significant'))])            
%             text(2.6,yl(2)*.9,[num2str(round(mean(mean(mean(dat(:,:,:,1) | isDown,3),2)),2,'significant'))])  
            ylabel('Proportion of neurons')            
            helper_saveandclosefig([thisdir 'PieEquivilant_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} '_1'])
        end


        if sum(ismember(figs_to_run,0))>0
            figure; hold on
            dat = zeros(size(isUp,1),size(isUp,2),size(isUp,3),3);
            dat(:,:,:,1) = islick;
            dat(:,:,:,2) = isUp; 
            dat(:,:,:,3) = isDown;     
            dat(repmat(sp.dem_num,[1 1 1 3])) = NaN;
            for imouse = 1:size(isUp,1)                
                plot(2+.1*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,1),3),'ok') 
                plot(3+.1*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,2),3),'ok')                 
                plot(4+.1*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,3),3),'ok')                 
            end
            yl = get(gca,'ylim');
            text(2.5,yl(2)*.99,[num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,1),3),2)),2,'significant'))])           
            text(3.5,yl(2)*.99,[num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,2),3),2)),2,'significant'))])                         
            text(4.5,yl(2)*.99,[num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,3),3),2)),2,'significant'))])   
            errorbar(2.5,nanmean(nanmean(nanmean(dat(:,:,:,1),3),2)),std(nanmean(nanmean(dat(:,:,:,1),3),2),1)./sqrt(size(dat(:,:,:,1),1)),'k','LineWidth',2)
            errorbar(3.5,nanmean(nanmean(nanmean(dat(:,:,:,2),3),2)),std(nanmean(nanmean(dat(:,:,:,2),3),2),1)./sqrt(size(dat(:,:,:,2),1)),'k','LineWidth',2)            
            errorbar(4.5,nanmean(nanmean(nanmean(dat(:,:,:,3),3),2)),std(nanmean(nanmean(dat(:,:,:,3),3),2),1)./sqrt(size(dat(:,:,:,3),1)),'k','LineWidth',2)            
            set(gca,'xtick',[2.25:4.25],'xticklabel',{'Lick-on';'Odor-On up';'Odor-On down'})            
            ylabel('Proportion of neurons')            
            helper_saveandclosefig([thisdir 'PieEquivilant_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} '_2'])
        end


        if sum(ismember(figs_to_run,[0 14]))>0
             figure; hold on
            dat = zeros(size(isUp,1),size(isUp,2),size(isUp,3),4);
            dat(:,:,:,1) = islickU;
            dat(:,:,:,2) = islickD;
            dat(:,:,:,3) = isUp; 
            dat(:,:,:,4) = isDown;     
            dat(repmat(sp.dem_num,[1 1 1 4])) = NaN;
            for imouse = 1:size(isUp,1)                
                for iii = 1:size(dat,4)
                    plot((iii+1)+.1*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,iii),3),'ok') 
                end
%                 plot(3+.1*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,2),3),'ok')                 
%                 plot(4+.1*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,3),3),'ok')                 
            end
            yl = get(gca,'ylim');
            for iii = 1:size(dat,4)
                text(iii+1.5,yl(2)*.99,[num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,iii),3),2)),2,'significant'))])                      
                errorbar(iii+1.5,nanmean(nanmean(nanmean(dat(:,:,:,iii),3),2)),std(nanmean(nanmean(dat(:,:,:,iii),3),2),1)./sqrt(size(dat(:,:,:,iii),1)),'k','LineWidth',2)
            end
%             text(2.5,yl(2)*.99,[num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,1),3),2)),2,'significant'))])           
%             text(3.5,yl(2)*.99,[num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,2),3),2)),2,'significant'))])                         
%             text(4.5,yl(2)*.99,[num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,3),3),2)),2,'significant'))])   
%             errorbar(2.5,nanmean(nanmean(nanmean(dat(:,:,:,1),3),2)),std(nanmean(nanmean(dat(:,:,:,1),3),2),1)./sqrt(size(dat(:,:,:,1),1)),'k','LineWidth',2)
%             errorbar(3.5,nanmean(nanmean(nanmean(dat(:,:,:,2),3),2)),std(nanmean(nanmean(dat(:,:,:,2),3),2),1)./sqrt(size(dat(:,:,:,2),1)),'k','LineWidth',2)            
%             errorbar(4.5,nanmean(nanmean(nanmean(dat(:,:,:,3),3),2)),std(nanmean(nanmean(dat(:,:,:,3),3),2),1)./sqrt(size(dat(:,:,:,3),1)),'k','LineWidth',2)            
            set(gca,'xtick',[2.25:5.25],'xticklabel',{'Lick-on up';'Lick-on down';'Odor-On up';'Odor-On down'})            
            ylabel('Proportion of neurons')            
            helper_saveandclosefig([thisdir 'PieEquivilant_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} '_2a'])
        end


        if sum(ismember(figs_to_run,0))>0
            figure; hold on
            dat = zeros(size(isUp,1),size(isUp,2),size(isUp,3),4);
            dat(:,:,:,1) = ~isDown & ~isUp & isOdorOff;
            dat(:,:,:,2) = ~isDown & isUp & isOdorOff; 
            dat(:,:,:,3) = ~isDown & isUp & ~isOdorOff;     
            dat(:,:,:,4) = isDown & ~isUp & ~isOdorOff;      
            dat(repmat(sp.dem_num,[1 1 1 4])) = NaN;
            for imouse = 1:size(isUp,1)                
                for iplot = 1:4
                    plot(iplot+1+.2*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,iplot),3),'ok') 
                end             
            end                       
            for iplot = 1:4
                errorbar(iplot+1.5,nanmean(nanmean(nanmean(dat(:,:,:,iplot),3),2)),std(nanmean(nanmean(dat(:,:,:,iplot),3),2),1)./sqrt(size(dat,1)),'k','LineWidth',2)
            end                    
            set(gca,'xtick',[2.5:5.5],'xticklabel',...
                {'Odor-off only';'Odor-off and odor-on';'Odor-on up only';'Odor-on down only'})
            yl = get(gca,'ylim');
            for iplot = 1:4
                text(iplot+1.5,yl(2)*.9,num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,iplot),3),2)),2,'significant')))
            end
            ylabel('Proportion of neurons')
            set(gcf,'Position',[1932         553        1040         498])      
            helper_saveandclosefig([thisdir 'PieEquivilant_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} '_3'])
        end


        if sum(ismember(figs_to_run,[0 15]))>0
            figure; hold on
             dat = zeros(size(isUp,1),size(isUp,2),size(isUp,3),5);
            dat(:,:,:,1) = ~isDown & ~isUp & isOdorOff;
            dat(:,:,:,2) = ~isDown & isUp & isOdorOff;
            dat(:,:,:,3) = isDown & ~isUp & isOdorOff;
            dat(:,:,:,4) = ~isDown & isUp & ~isOdorOff;
            dat(:,:,:,5) = isDown & ~isUp & ~isOdorOff;              
            dat(repmat(sp.dem_num,[1 1 1 5])) = NaN;
            for imouse = 1:size(isUp,1)                
                for iplot = 1:5
                    plot(iplot+1+.2*imouse+rand(size(isUp,2),1)/20,nanmean(dat(imouse,:,:,iplot),3),'ok') 
                end   
            end
            for iplot = 1:5
                errorbar(iplot+1.5,nanmean(nanmean(nanmean(dat(:,:,:,iplot),3),2)),std(nanmean(nanmean(dat(:,:,:,iplot),3),2),1)./sqrt(size(dat,1)),'k','LineWidth',2)
            end            
            yl = get(gca,'ylim');
            for iplot = 1:5
                text(iplot+1.5,yl(2)*.9,num2str(round(nanmean(nanmean(nanmean(dat(:,:,:,iplot),3),2)),2,'significant')))
            end
            set(gca,'xtick',[2.5:6.5],'xticklabel',...
                {'Odor-off only';'Odor-off and odor-on up';'Odor-off and odor-on down';'Odor-on up only';'Odor-on down only'})
            ylabel('Proportion of neurons')
            set(gcf,'Position',[1932         553        1040         498])            
            helper_saveandclosefig([thisdir 'PieEquivilant_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana} '_4'])
        end
    end        
end