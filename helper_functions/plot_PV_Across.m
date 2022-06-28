function plot_PV_Across(corr_repeat_trial_session,all_labels,dirs,itype,smd)
thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'PV_Across/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

if itype==2
    itt2 = [1 3 4; 1 2 4; 1 2 3];
elseif itype==1
    itt2 = repmat(1:4,[3 1]); % to look at how probes deal with probes
end

for itarget = length(all_labels.labss):-1:1
    for iplottype = 2 %1:2 
      for iraw = 1:2
          for icor = 1:3
              if iplottype == 1
                  plotind = 6; plotlabel = 'Trial'; plotlabeloth = 'Session'; plotind2 = 7; binsz = 10;
              elseif iplottype == 2
                  plotind = 7; plotlabel = 'Session'; plotlabeloth = 'Trial'; plotind2 = 6; binsz = 1;
              end

             OGdat = corr_repeat_trial_session;

             if itarget<3
                indall = OGdat(:,5)==itarget;
             else             
                indall = OGdat(:,5)>2;
             end

             if icor ==1
                 corlab = 'AllTrials';
             elseif icor == 2
                 corlab = 'CorrectOnly';
                 indall = indall & OGdat(:,9)==1;
             elseif icor == 3             
                 corlab = 'IncorrectOnly';
                 indall = indall & OGdat(:,9)==0;
             end

             dat = OGdat(indall,:);
                
             clear ff
             ffylim = NaN(2,size(itt2,2));
             figure; hold on
             for itt = 1:size(itt2,2)
                if iraw == 1
                    datuse = dat(:,itt2(itarget,itt)); rawlab = 'Raw';
                elseif iraw==2
                       [~,~,stats] = anovan(dat(:,itt2(itarget,itt)),dat(:,[plotind2 8:9]),'varnames',{plotlabeloth,'Mouse','Rewarded'},'continuous',1,'display','off');    
                        datuse = stats.resid; rawlab = 'Residuals';
                        if size(datuse,1)==1
                            datuse = datuse';
                        end
                        if iplottype == 2
                          [~,tbl] = anovan(dat(:,itt2(itarget,itt)),dat(:,[plotind plotind2 8:9]),'varnames',{plotlabel,plotlabeloth,'Mouse','Rewarded'},'continuous',[1 2],'display','off');    
                           writecell(tbl,[thisdir all_labels.typelab{itype} '_PVPlot_' all_labels.labss{itarget} '_Over' plotlabel '_' rawlab '_' corlab all_labels.addon '_ANOVA_' all_labels.ittlabs{itt2(itarget,itt)} '.txt'],'Delimiter','tab')    
                        end
                end
                ff(itt) = subplot(1,size(itt2,2),itt); hold on;
                    if iplottype==1
                        plot(dat(:,plotind),datuse,'o','Color',[.5 .5 .5])
                    elseif iplottype==2
                        plot(dat(:,plotind),datuse,'.k')
%                             violinplot1(datuse,dat(:,plotind),'EdgeColor',[.5 .5 .5],'BoxColor',[.5 .5 .5],'ViolinColor',[.5 .5 .5],'MedianColor',[.3 .3 .3]);
                    end
                dat1 = ceil((dat(:,plotind)-min(dat(:,plotind))+.001)/binsz);
                if sum(diff(sort(unique(dat1)))>1)==0
                    
                       
                    m1 = splitapply(@nanmean,datuse,dat1);                        
                    m2 = smoothdata([m1(end:-1:1) m1 m1(end:-1:1)],'gaussian',smd);
                    m = m2(length(m1)+1:length(m1)*2)';
                    datss = splitapply(@nanstd,datuse,dat1);
                    h = hist(dat1,1:length(m));
                    sem1 = datss./sqrt(h');                        
                    sem2 = smoothdata([sem1(end:-1:1) sem1 sem1(end:-1:1)],'gaussian',smd);
                    sem = sem2(length(sem1)+1:length(sem1)*2)';
                    rev = (m+sem)'; fwd = (m-sem)';
                    tbs = (1:length(m))*binsz;
                    p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
                    p1.FaceAlpha=.2; p1.EdgeAlpha=0;   % p1.FaceColor = [0 0 1];
                    plot(tbs,m,'k');
                    
                else
    %                         if iplottype==2
    %                              plot(dat(:,plotind),datuse,'o','Color',[.5 .5 .5])
    %                         end
                end
                 xlim([1 max(dat(:,plotind))])
                [r,p] = corr(dat(:,plotind),datuse,'rows','complete');    
                 numshuffhere = 1000;
                 r_shuff = NaN(numshuffhere,1);
                 datjnk = dat(:,plotind);
                 for ishuff = 1:numshuffhere                             
                    r_shuff(ishuff) = corr(datjnk(randperm(length(datjnk))),datuse,'rows','complete');    
                 end
                 p2 = (1+sum(abs(r_shuff)>=abs(r)))./(1+all_labels.numshuff);
                     
                N = sum(~isnan(datuse) & ~isnan(dat(:,plotind)));
                tt = title({[all_labels.ittlabs{itt2(itarget,itt)} ', N = ' num2str(N)];['R = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant')) ', ShuffP = ' num2str(round(p2,2,'significant'))]});
                if p<.05 && r>0
                    tt.Color = 'r';
                elseif p<.05 && r<0
                    tt.Color = 'b';
                end
                 xlabel(plotlabel)
                 ylabel(all_labels.ittlabs{itt2(itarget,itt)})
                ffylim(:,itt) = get(gca,'ylim');
                xlimy = get(gca,'xlim');
                set(gca,'xlim',[0 xlimy(2)+1])
             end
             suptitle({[all_labels.typelab{itype} ' experiment, ' all_labels.labss{itarget} ' Trials'];['Over ' plotlabel 's, ' rawlab ', ' corlab]})
             set(gcf,'Position',[ 1959         173        1010         650])
             for itt = 1:size(itt2,2)
                 set(ff(itt),'ylim',[min(ffylim(1,:)) max(ffylim(2,:))])
             end
             helper_saveandclosefig([thisdir all_labels.typelab{itype} '_PVPlot_' all_labels.labss{itarget} '_Over' plotlabel '_' rawlab '_' corlab all_labels.addon])
             
             
          end
      end
    end
end
  