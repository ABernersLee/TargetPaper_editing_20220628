function helper_plot_lick_lat(lickall,dirs)
    %%

%from start of decision period
lick_cutoff_ms = 500;


lickall2 = lickall(:,1); 
lickall2(lickall(:,1)>100) = 3; %probe
lickall2(lickall(:,1)>2 & lickall(:,1)<100) = 4; %NT repeat

if ~isfolder([dirs.figdir '\Behavior\NeuralDataSessions\Licks\'])
    mkdir([dirs.figdir '\Behavior\NeuralDataSessions\Licks\'])
end

% lickall = trialtype, latency of first lick, mouse, accuracy, session, trial #, type of session   

lickacc = cat(2,lickall(:,3),lickall(:,5),lickall(:,4),lickall(:,2),lickall2(:,1),lickall(:,6:7));
% lickacc(lickacc(:,7)~=1,:) = []; %probe sessions only
lickacc(lickacc(:,4)<=lick_cutoff_ms,:) = [];
tlabs = {'Target';'NT';'Other'};


%% by trial type

clump = 8; %12 %venki asked to clump together the later latencies

clabs = {'Incorrect','Correct'};
binsize = 100; stnd = [500 2200];
ind = stnd(1):binsize:stnd(2);
% ind = 400:binsize:stnd(2);
Ms = NaN(length(ind),4); SD = Ms; hs = NaN(length(ind)-1,4);
for itt = 1:3
    lickacc2 = lickacc(lickacc(:,5)==itt,:);
    [h,~,c] = histcounts(lickacc2(:,4),ind);
    c(c>(clump-1)) = clump;
    hs(:,itt) = h;
    M = splitapply(@mean,lickacc2(c>0,3),c(c>0));
    Ms(unique(c(c>0)),itt) = M;
    S = splitapply(@std,lickacc2(c>0,3),c(c>0));
    D = splitapply(@sum,ones(sum(c>0),1),c(c>0));
    SD(unique(c(c>0)),itt) = S./sqrt(D);
end


% [h,~] = histc(lickacc(:,4),ind);
vc = [0 .6 0;0 0 .6; .6 0 0;.6 0 .6];
figure; hold on;
set(gcf,'Position',[87         135        1621         786])

% subplot(2,4,[3:4 7:8]); hold on
subplot(2,3,[2:3 5:6]); hold on
yyaxis left
% bar(h./sum(h),'FaceColor',[.5 .5 .5],'EdgeColor','w')
for itt = 1:3
    bar(hs(:,itt)./sum(hs(:,itt)),'FaceColor',vc(itt,:),'EdgeColor','w','FaceAlpha',.3)
end
yl = get(gca,'ylim');
set(gca,'ylim',[yl(1) yl(2)+(range(yl)*.5)])
ylabel('Probability')
yyaxis right
for itt = 1:3
%     errorbar(Ms(:,itt),SD(:,itt),'Color',vc(itt,:),'LineStyle','-')
    tbs = 1:length(ind);
%     tbs = [1:length(ind)-1 (length(ind)*2+clump)/2];
    rev = (Ms(:,itt)+SD(:,itt))'; fwd = (Ms(:,itt)-SD(:,itt))';
    tbs(isnan(rev)) = []; fwd(isnan(rev)) = []; rev(isnan(rev)) = [];  
    tbs = [tbs(1:end-1)  (tbs(end)+length(ind)-1)/2]; %((length(ind)-1)-clump)/2+clump];
    p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
    p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = vc(itt,:);
    plot(tbs,Ms(~isnan(Ms(:,itt)),itt),'.-','Color',vc(itt,:),'MarkerSize',20)
end
ylabel('Proportion correct')
set(gca,'xtick',1:5:length(ind),'xticklabel',ind(1:5:length(ind)))
xlabel('Lick latency (ms)')
legend('Target','NT','Probe','Location','Best')


for it = 1:2
%     subplot(2,4,1+(it-1)*4); hold on
    subplot(2,3,1+(it-1)*3); hold on
    ind = lickacc(:,5)==it;
    histogram(lickacc(lickacc(:,3)==1 & ind,4),'FaceColor',vc(it,:),'EdgeColor','w','FaceAlpha',.3)
    histogram(lickacc(lickacc(:,3)==0 & ind,4),'FaceColor','k','EdgeColor','w','FaceAlpha',.3)
    xlabel('Latency (ms)')
    ylabel('Count')    
    yl2 = get(gca,'ylim');
    plot([median(lickacc(lickacc(:,3)==1 & ind,4),'omitnan') median(lickacc(lickacc(:,3)==1 & ind,4),'omitnan')],yl2,'-','color',vc(it,:))
    plot([median(lickacc(lickacc(:,3)==0 & ind,4),'omitnan') median(lickacc(lickacc(:,3)==0 & ind,4),'omitnan')],yl2,'-','color','k')
%     text(1100,yl2(2)-range(yl2)/2,['CorrectN = ' num2str(sum(lickacc(:,3)==1 & ind)) ', InCorrN = ' num2str(sum(lickacc(:,3)==0 & ind))])
    legend('Correct','Incorrect','Location','Best')
    title(tlabs{it})
    
%     [p,tbl,stats] = anova1(lickacc(ind,4),lickacc(ind,3),'off');
    p = anovan(lickacc(ind,4),lickacc(ind,[3 1:2 6]),'Display','off','varnames',{'Accuracy','Mouse','Day','Trial#'},'continuous',3:4);
    
    text(1500,yl2(2)-range(yl2)/2,['p = ' num2str(round(p(1),2,'significant'))])
end

% for it = 1:2
%     subplot(2,4,2+(it-1)*4); hold on
%     ind = lickacc(:,3)==(it-1);
%     for itt = 1:3
%         histogram(lickacc(lickacc(:,5)==itt & ind,4),'FaceColor',vc(itt,:),'EdgeColor','w','FaceAlpha',.3)
%     end
%     xlabel('Latency (ms)')
%     ylabel('Count')    
%     yl2 = get(gca,'ylim');
%     for itt = 1:3
%         plot([median(lickacc(lickacc(:,5)==itt & ind,4),'omitnan') median(lickacc(lickacc(:,5)==itt & ind,4),'omitnan')],yl2,'-','color',vc(itt,:))
%     end
% %     text(1100,yl2(2)-range(yl2)/2,['TargetN = ' num2str(sum(lickacc(:,5)==1 & ind)) ', NtN = ' num2str(sum(lickacc(:,5)==2 & ind))])
%     legend('Target','NT','Probe','Location','Best')
%     title(clabs{it})
%     
%     ind = ind & lickacc(:,5)<4;
% %     [p,tbl,stats] = anova1(lickacc(ind,4),lickacc(ind,5),'off');
%     [p,~,stats] = anovan(lickacc(ind,4),lickacc(ind,[5 1:2 6]),'Display','off','varnames',{'Trialtype','Mouse','Day','Trial#'},'continuous',3:4);
%     pp = multcompare(stats,'display','off','Dimension',1);    
%     text(1500,yl2(2)-range(yl2)/2,['p = ' num2str(round(p(1),2,'significant'))])
% end

%%
helper_saveandclosefig([dirs.figdir ...
    '\Behavior\NeuralDataSessions\Licks\SpeedAccuracyTradeoff'])
 