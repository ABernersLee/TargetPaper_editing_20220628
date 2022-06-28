function helper_plot_lick(lickall,mice,mousenames,dirs,figs_to_run)
    
labcol = [0 .9 0;0 0 .9];

%from start of decision period
lick_cutoff_ms = 500;
lickall(lickall(:,2)<=lick_cutoff_ms,:) = [];
% lickall(:,2) = lickall(:,2)-lick_cutoff_ms;

Ls = [111 222 333; 11 22 33];
lickall2 = lickall(:,1); 
lickall2(lickall(:,1)>2 & lickall(:,1)<100) = 3;
lickall2(lickall(:,1)>100) = 4;

if ~isfolder([dirs.figdir '\Behavior\NeuralDataSessions\Licks\'])
    mkdir([dirs.figdir '\Behavior\NeuralDataSessions\Licks\'])
end

% trialtype, latency of first lick, mouse, accuracy, session, trial #, type of session 

% figure 1 - lick latencies
ind = lickall(:,4)==1;
[~,~,stats1] = anovan(lickall(ind,2),...
    [lickall2(ind) lickall(ind,[3 5 6])],'display','off',...
    'varnames',{'Trialtype';'Mouse';'Session';'Trial'},'continuous',3:4);
[stats2,s3] = multcompare(stats1,'Dimension',1,...
    'ctype','bonferroni','display','off');

if sum(ismember(figs_to_run,[0 17]))>0
    figure; hold on

    errorbar(1:size(s3,1),s3(:,1),s3(:,2),'CapSize',0,...
        'Color','k','Marker','o','LineStyle','none','MarkerSize',10);
    xlim([.5 size(s3,1)+.5])
    ylabel('Marginal mean lick latency')
    yy = get(gca,'ylim');
    for istat = 1:size(stats2,1)
        if stats2(istat,6)<0.05
    
            plot([stats2(istat,1) stats2(istat,2)], ...
                [yy(2)+istat*(range(yy)/15) yy(2)+istat*(range(yy)/15)],'-k')  
    
            if stats2(istat,6)<0.01 &&  stats2(istat,6)>=0.001
               text(mean(stats2(istat,1:2)), ...
                   yy(2)+istat*(range(yy)/15)+range(yy)/75,'**','color','k')
            elseif stats2(istat,6)<0.001
               text(mean(stats2(istat,1:2)), ...
                   yy(2)+istat*(range(yy)/15)+range(yy)/75,'***','color','k')
            else
               text(mean(stats2(istat,1:2)), ...
                   yy(2)+istat*(range(yy)/15)+range(yy)/75,'*','color','k')
            end
    
        end
    end
    set(gca,'xtick',1:size(s3,1), ...
        'xticklabel',{'Target';'Non-target';'NT repeat';'Probe'})
    set(gcf,'Position',[2193         385         243         415])
    helper_saveandclosefig([dirs.figdir ...
        '\Behavior\NeuralDataSessions\Licks\AcrossAllSessionsUsing_Lick_ANOVA_latency'])
end

%%% Figure 5
if sum(ismember(figs_to_run,[0 5]))>0
    itype = 2;
    pcol1 = [.7 1 .7; .7 .7 1; .5 .5 .5];
    pcol = [0 .9 0; 0 0 .6;0 0 0];
    mdat2 = [];
    clear p
    figure; hold on;
    for itrialtype = 1:3
        subplot(3,1,itrialtype); hold on
        mdat = NaN(30,length(mice));
        for imouse = 1:length(mice)
            ind = lickall(:,3)==mice(imouse) & lickall(:,7)==itype & lickall2(:,1)==itrialtype & lickall(:,4)==1;     
    %         newsessions = unique(lickall(ind,5));
            for isession = 1:max(lickall(ind,5))
                mdat(isession,imouse) = mean(lickall(ind & lickall(:,5)==isession,2),'omitnan');
            end
            plot(unique(lickall(ind,5)),mdat(~isnan(mdat(:,imouse)),imouse),'o-','Color',pcol1(itrialtype,:))
        end
       
        p(itrialtype) = errorbar(mean(mdat,2,'omitnan'), ...
            std(mdat,[],2,'omitnan')./sqrt(sum(~isnan(mdat),2,'omitnan')), ...
            '-','Color',pcol(itrialtype,:)); 
    
        axis tight
        xl = get(gca,'xlim');
        xlim([xl(1)-1 xl(2)+1])
        ylim([-500 2000])
        set(gca,'xtick',5:5:25,'ytick',0:500:1500)
        xlabel('Session')
        ylabel('Lick latency (ms)')
        mdat2 = cat(3,mdat2,mdat);
    %     legend(p,'target','nontarget','nt repeats','Location','best')
    end
    mdat3 = mean(mdat2,3,'omitnan');
    B1 = (ones(size(mdat3,2),1)*(1:size(mdat3,1)))';
    A1 = ones(size(mdat2,1),1)*(1:size(mdat2,2));
    B = repmat(B1,[1 1 3]);
    A = repmat(A1,[1 1 3]);
    C1 = ones(size(mdat2,1),size(mdat2,2));
    C = cat(3,C1,C1*2,C1*3);
    N = sum(sum(~isnan(mdat3)));
    [ps,tbl,~] = anovan(mdat2(:),...
        [A(:) B(:) C(:)],'display','off',...
        'varnames',{'Mouse';'Session';'TrialType'},'continuous',2);
    [r,p] = corr(mdat3(:),B1(:),'rows','complete');
    set(gcf,'units','centimeters'); set(gcf,'Position',[0 0 9 24])     
    trp = helper_text_rp(r,p,2);
    suptitle({['N = ' num2str(N) ', NT repeats'];[trp ' ANOVAp = ' num2str(ps(2))]})
    helper_saveandclosefig([dirs.figdir ...
        '\Behavior\NeuralDataSessions\Licks\AcrossSessions_Latency'])
    writecell(tbl,[dirs.figdir ...
        '\Behavior\NeuralDataSessions\Licks\AcrossSessions_Latency.xlsx'])
end

tt_labs = {'T';'NT';'NT-R';'P'};

if sum(ismember(figs_to_run,[0 17 18]))>0
    
     if sum(ismember(figs_to_run,[0 17]))>0
        % supplemental figure 7
        fig_a = figure; hold on; 
        clear aa; yys = NaN(length(mice),2);
     end

    for imouse = 1:length(mice)

        if sum(ismember(figs_to_run,[0 17]))>0
            figure(fig_a)
            aa(imouse) = subplot(1,4,imouse); hold on
        end
            
        ind = lickall(:,4)==1 & lickall(:,3)==imouse;
        [~,~,stats1] = anovan(lickall(ind,2),[lickall2(ind) lickall(ind,5)],...
            'display','off','varnames',{'Trialtype';'Session'},'continuous',2);
        [stats2,s3] = multcompare(stats1,'Dimension',1,'ctype','bonferroni',...
            'display','off');
    
        ts = unique(lickall2(ind & ~isnan(lickall2)));
        ll2 = unique(lickall(ind,5));
        A = ceil(sqrt(length(ll2)));
        B = ceil(length(ll2)./A);
        
        if (sum(ismember(figs_to_run,18))>0 && strcmp('Q',mousenames{imouse})) || sum(ismember(figs_to_run,0))>0
            % supplemental figure 8 (mouse Q)
            fig_b = figure; hold on 
            for isess = 1:length(ll2)        
                subplot(A,B,isess); hold on
                for itt = 1:length(ts)
                    thisind = lickall(:,4)==1 & lickall(:,3)==imouse & lickall(:,5)==ll2(isess) & lickall2==ts(itt);
                    if sum(thisind)>0
                        Violin(lickall(thisind,2),ts(itt),'ViolinColor',[.5 .5 .5]);
        %                     errorbar(ts(itt),mean(lickall(thisind,2)),std(lickall(thisind,2))./sqrt(sum(thisind)),'k')
                %         Violin(lickall(lickall2==ts(itt),2),itt,'ViolinColor',[.5 .5 .5]);
                    end
                end
                xlim([0.5 4.5])
                title(['Session ' num2str(ll2(isess))])
                set(gca,'xtick',ts,'xticklabel',tt_labs(ts))
                ylabel('Lick latency (ms)')
            end
            set(fig_b,'Position',[679          75        1082         883])
            helper_saveandclosefig([dirs.figdir ...
                '\Behavior\NeuralDataSessions\Licks\AllSessions_mouse' mousenames{imouse}])
        end
    
    
        if sum(ismember(figs_to_run,[0 17]))>0

            figure(fig_a)
            errorbar(1:size(s3,1),s3(:,1),s3(:,2),...
                'CapSize',0,'Color','k','Marker','.','LineStyle','none'...
                ,'MarkerSize',10);
        
            xlim([.5 size(s3,1)+.5])
            ylabel('Marginal mean lick latency (ms)')
            set(gca,'xtick',1:size(s3,1),'xticklabel',tt_labs)
        
            title(['Mouse ' mousenames{imouse}])
            ylabel('Cumulative probability')
            xlabel('Lick latency (ms)')
            set(gca,'FontSize',12)
            yy = get(gca,'ylim');    
            mu = (range(yy)/24); %12);
            for istat = 1:size(stats2,1)
                if stats2(istat,6)<0.05
                    plot([stats2(istat,1) stats2(istat,2)],...
                        [yy(2)+istat*mu yy(2)+istat*mu],'-k')        
                    if stats2(istat,6)<0.01 &&  stats2(istat,6)>=0.001
                       text(mean(stats2(istat,1:2)),yy(2)+istat*mu+2,...
                           '**','color','k')
                    elseif stats2(istat,6)<0.001
                       text(mean(stats2(istat,1:2)),yy(2)+istat*mu+2,...
                           '***','color','k')
                    else
                       text(mean(stats2(istat,1:2)),yy(2)+istat*mu+2, ...
                           '*','color','k')
                    end
                end
            end
%             axis tight
%             set(gca,'xlim',[.5 size(s3,1)+.5])
            yys(imouse,:) = get(gca,'ylim');
        end
    end

    if sum(ismember(figs_to_run,[0 17]))>0
        %to eq. axes
%         for imouse = 1:length(mice)
%             set(aa(imouse),'ylim',[min(yys(:,1)) max(yys(:,2))])
%         end
%         set(gcf,'Position',[ 205         686        1318         216])
        set(gcf,"Position",[205         290        1318         612])
        helper_saveandclosefig([dirs.figdir ...
            '\Behavior\NeuralDataSessions\Licks\AcrossAllSessionsUsing_Lick_ANOVA_mice'])
    
    
        % supplemental figure 7
        figure; hold on; 
        px = NaN(length(mice),2);
        clear aa
        for imouse = 1:length(mice)
            aa(imouse) = subplot(1,4,imouse); hold on
            clear p
            for itt = 1:2
                ind = lickall(:,3)==imouse & lickall(:,1)==itt;                    
                p(itt) = cdfplot(lickall(ind,2));
                p(itt).Color = labcol(itt,:);
                p(itt).LineWidth = 2;
            end     
        
            for itype = 1:2                                          
                numtt = 5-(itype-1);
                for itt = 3:numtt
                    ind = lickall(:,3)==imouse & lickall(:,1)==Ls(itype,itt-2);
        
                    if sum(ind)>0
                        p(itype+2) = cdfplot(lickall(ind,2));
                        p(itype+2).Color = labcol(2,:);
                        if itype == 1
                            p(itype+2).LineStyle = '--';
                        end
                    end
                end
            end
            title(['Mouse ' mousenames{imouse}])
            ylabel('Cumulative probability')
            xlabel('Lick latency (ms)')
            set(gca,'FontSize',12)
            px(imouse,:) = get(gca,'xlim');
            if imouse == 2
                legend(p,{'Target';'NT';'Probes';'Repeats'},'Location','Best')
            end
        end
        for imouse = 1:length(mice)
            aa(imouse).XLim = [min(px(:,1)) max(px(:,2))];
        end

        set(gcf,'Position',[ 205         686        1318         216])
        helper_saveandclosefig([dirs.figdir ...
            '\Behavior\NeuralDataSessions\Licks\AcrossAllSessionsUsing_Lick'])
    end
end
