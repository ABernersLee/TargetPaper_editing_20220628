function plot_Selectivity_Comparisons(singledat,selectivity_dat,selectivity_params,dirs,figs_to_run) 

% right now only one figure comes out of this for figure 4

od = singledat.toplotOdorLick(:,1:6000);

sp = selectivity_params;
sd = selectivity_dat;

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'Comparisons/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end
     
  %%%%% Comparisons between categories of neurons, between timing too
times = sp.odorwindowlabel(:,1)>=0 & sp.odorwindowlabel(:,2)<=2000;

for iana = [1 3]
   
    is_target = sd.all_Neurons(:,5);                
    is_nontarget = sd.all_Neurons(:,6);    
    isUp = sd.all_Neurons(:,3);            
    isDown = sd.all_Neurons(:,4);
    isOdorOff = sd.all_Neurons(:,13);
    is_target_off = sd.all_Neurons(:,15);
    is_nontarget_off = sd.all_Neurons(:,16);
    isSigOdor = sd.all_Neurons(:,1);

%     if iana==2
%         isUp = isUp & sd.neuron_info(:,7)==1;                    
%         isDown = isDown & sd.neuron_info(:,7)==1;
%         is_target = is_target & sd.neuron_info(:,7)==1;
%         is_nontarget = is_nontarget & sd.neuron_info(:,7)==1;
%         isOdorOff = isOdorOff & sd.neuron_info(:,7)==1;
%         is_target_off = is_target_off & sd.neuron_info(:,7)==1;
%         is_nontarget_off = is_nontarget_off & sd.neuron_info(:,7)==1;
%         isSigOdor = isSigOdor & sd.neuron_info(:,7)==1;
%     end
     
            
    odorbins = squeeze(singledat.odorsel(:,times,1,:));
    upCount = sum(odorbins(:,:,1)>0 & odorbins(:,:,2)<(sp.p_cutoff/size(odorbins,2)),2);
%     downCount = sum(odorbins(:,:,1)<0 & odorbins(:,:,2)<(sp.p_cutoff/size(odorbins,2)),2);
    
    upHist = odorbins(:,:,1)>0 & odorbins(:,:,2)<(sp.p_cutoff/size(odorbins,2));
%         downHist = sum(odorbins(:,:,1)<0 & odorbins(:,:,2)<(sp.p_cutoff/size(odorbins,2)),1);
    
    targetbins = squeeze(singledat.odorsel(:,times,3,:));
    upCountT = sum(targetbins(:,:,2)<(sp.p_cutoff/size(targetbins,2)),2);     
    
    upCountTU = sum(targetbins(:,:,2)<(sp.p_cutoff/size(targetbins,2)) & odorbins(:,:,1)>0 & odorbins(:,:,2)<(sp.p_cutoff/size(odorbins,2)),2);     
    upHistT = targetbins(:,:,2)<(sp.p_cutoff/size(targetbins,2));
            
    % Chi-squared test between up and down modulated and modulated by target at all
    % Observed data
    n1 = sum((is_target | is_nontarget) & isDown); N1 = sum(isDown);
    n2 = sum((is_target | is_nontarget) & isUp); N2 = sum(isUp);
    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2);
    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    % Chi-square test, by hand
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    [~,pp0,PPstats] = chi2gof([1 2 3 4],'freq',observed,'expected',expected,'ctrs',[1 2 3 4],'nparams',2);
    PPstats.Pvalue = pp0;

    D = [sum(~is_target & ~is_nontarget & isDown) sum(is_target & isDown) sum(is_nontarget & isDown)]./(sum(isDown));
    U = [sum(~is_target & ~is_nontarget & isUp) sum(is_target & isUp) sum(is_nontarget & isUp)]./(sum(isUp));
    Ubreif = [sum(~is_target & ~is_nontarget & isUp & upCount<sp.medcut ) sum(is_target & isUp & upCount<sp.medcut ) sum(is_nontarget & isUp & upCount<sp.medcut )]./(sum(isUp & upCount<sp.medcut ));
    Ulong = [sum(~is_target & ~is_nontarget & isUp & upCount>=sp.medcut ) sum(is_target & isUp & upCount>=sp.medcut ) sum(is_nontarget & isUp & upCount>=sp.medcut )]./(sum(isUp & upCount>=sp.medcut));
    O = [sum(~is_target_off & ~is_nontarget_off & isOdorOff) sum(is_target_off & isOdorOff) sum(is_nontarget_off & isOdorOff)]./(sum(isOdorOff));

    % Chi-squared test between long and briefly modulated and T.Nt preferring
    % Observed data
    n1 = sum(is_target & upCount<sp.medcut  ); N1 = sum(upCount<sp.medcut  & (is_target | is_nontarget) );
    n2 = sum(is_target & upCount>=sp.medcut ); N2 = sum(upCount>=sp.medcut  & (is_target | is_nontarget) );
    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2);
    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    % Chi-square test, by hand
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    [~,p00] = chi2gof([1 2 3 4],'freq',observed,'expected',expected,'ctrs',[1 2 3 4],'nparams',2);


    if (sum(ismember(figs_to_run,[0 4]))>0 && iana==1) || (sum(ismember(figs_to_run,0))>0)
        figure; hold on
%         bar([D; U; Ubreif; Ulong; O],'stacked')
        bar([D; U],'stacked')
        n1 = sum((is_target | is_nontarget) & isUp);
        n2 = sum(isUp);
        nn1 = sum((is_target | is_nontarget) & isDown);
        nn2 = sum(isDown);
        p01 = 1-binocdf(n1-1,n2,sp.p_cutoff*2);
        p02 = 1-binocdf(nn1-1,nn2,sp.p_cutoff*2);
        set(gca,'xtick',1:2,'xticklabel',{'Down';'Up'})
        text(2.5,.2,{'Chi-squared test of Up vs Down'; 'being modulated by target at all';['p = ' num2str(round(pp0,2,'significant'))]})
        text(2.5,.85,{'Binomial probability of Up';'cells being modulated,';[num2str(n1) ' of ' num2str(n2) ', p = ' num2str(round(p01,2,'significant'))]})
        text(2.5,.55,{'Binomial probability of Down';'cells being modulated,';[num2str(nn1) ' of ' num2str(nn2) ', p = ' num2str(round(p02,2,'significant'))]})        
        xlim([0 7])
        legend({'Not mod.','Target','NonTarget'},'Location','northeast')
        ylabel('Proportion of cells')        
        helper_saveandclosefig([thisdir 'PropSigUpDownTNT_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_' sp.ana_labs{iana}])
        
        PP(1,:) = fieldnames(PPstats);
        PP(2,:) = struct2cell(PPstats);
        writecell(PP,[thisdir 'PropSigUpDownTNT_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_' sp.ana_labs{iana} '_stats.xlsx'])
    end


    if sum(ismember(figs_to_run,0))>0
        figure; hold on
        bar([Ubreif; Ulong],'stacked')        
        set(gca,'xtick',1:2,'xticklabel',{'Up brief';'Up long'})
        text(2.5,.2,{'Chi-squared test';'of Brief vs Long';['p = ' num2str(round(p00,2,'significant'))]})
        xlim([0 4])
        legend({'Not mod.','Target','NonTarget'},'Location','northeast')
        ylabel('Proportion of cells')

        
        helper_saveandclosefig([thisdir 'PropSigBriefLongTNT_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_medcut' num2str(sp.medcut) '_' sp.ana_labs{iana}])
        
        figure; hold on
        subplot(5,2,3); hold on
        histogram(upCount(isUp & ~is_target & ~is_nontarget & upCount>0),'Normalization','probability')
        histogram(upCount(isUp & (is_target | is_nontarget) & upCount>0),'Normalization','probability')
        p1 = ranksum(upCount(isUp & ~is_target & ~is_nontarget & upCount>0),upCount(isUp & (is_nontarget | is_target) & upCount>0));
        ylabel('Proportion of cells')
        xlabel('Number of bins sig up')
        legend({'Not T/NT modulated';'T/NT Modulated'},'Location','best')
        title(['Up, p = ' num2str(round(p1,2,'significant'))])

        subplot(5,2,4); hold on
        histogram(upCount(isUp & is_target & upCount>0))
        histogram(upCount(isUp & is_nontarget & upCount>0))
        p2 = ranksum(upCount(isUp & is_nontarget & upCount>0),upCount(isUp & is_target & upCount>0));
        ylabel('Proportion of cells')
        xlabel('Number of bins sig up')
        legend({'Target preference'; 'Non-target preference'},'Location','best')
        title(['Up, p = ' num2str(round(p2,2,'significant'))])                       

        subplot(5,2,5); hold on
        dat1 = [upCountT(isUp & is_target & upCountT>0)];
        dat2 = [upCountT(isUp & is_nontarget & upCountT>0)];
        histogram(dat1,'Normalization','probability')
        histogram(dat2,'Normalization','probability')
        p2 = ranksum(dat1,dat2);
        ylabel('Proportion of cells')
        xlabel('Number of bins sig Target/NT')
        legend({'Target preference'; 'Non-target preference'},'Location','best')
        title(['Up, p = ' num2str(round(p2,2,'significant'))])
        
         subplot(5,2,6); hold on
        dat1 = upCountTU(isSigOdor & is_target & upCountTU>0);
        dat2 = upCountTU(isSigOdor & is_nontarget & upCountTU>0);
        histogram(dat1,'Normalization','probability')
        histogram(dat2,'Normalization','probability')
        p1 = ranksum(dat1,dat2);
        ylabel('Proportion of cells')
        xlabel('Number of bins sig up and Target/NT')
        legend({'Target preference'; 'Non-target preference'},'Location','best')
        title(['Up, p = ' num2str(round(p1,2,'significant'))])
        
         subplot(5,2,7); hold on
        dat1 = upHist(isUp & ~is_target & ~is_nontarget & upCount>0,:);
        dat2 = upHist(isUp & (is_target | is_nontarget) & upCount>0,:);
        bar(mean(dat1,1),'FaceColor','none','EdgeColor','b')%,'Normalization','probability')
        bar(mean(dat2,1),'FaceColor','none','EdgeColor','k')%,'Normalization','probability')                
        p2 = anovan([dat1(:);dat2(:)],[[ones(size(dat1(:)));2*ones(size(dat2(:)))] [repmat([1:size(dat1,2)]',[size(dat1,1) 1]);repmat([1:size(dat2,2)]',[size(dat2,1) 1])]],'display','off','model','interaction');       
        ylabel('Proportion of cells')
        xlabel('Number of cells sig per bin')
        legend({'Not T/NT modulated';'T/NT Modulated'},'Location','best')
        title(['Up, p = ' num2str(round(p2(3),2,'significant'))])

        subplot(5,2,8); hold on
        dat1 = upHist(isUp & is_nontarget & upCount>0,:);
        dat2 = upHist(isUp & is_target & upCount>0,:);        
        bar(mean(dat1,1),'FaceColor','none','EdgeColor','b')%,'Normalization','probability')
        bar(mean(dat2,1),'FaceColor','none','EdgeColor','k')%,'Normalization','probability')                
        p2 = anovan([dat1(:);dat2(:)],[[ones(size(dat1(:)));2*ones(size(dat2(:)))] [repmat([1:size(dat1,2)]',[size(dat1,1) 1]);repmat([1:size(dat2,2)]',[size(dat2,1) 1])]],'display','off','model','interaction');      
        ylabel('Proportion of cells')
        xlabel('Number of bins sig up')
        legend({'Target preference'; 'Non-target preference'},'Location','best')
        title(['Up, p = ' num2str(round(p2(3),2,'significant'))])
                
        subplot(5,2,9); hold on
        dat1 = upHistT(isSigOdor & is_target & upCountT>0,:);
        dat2 = upHistT(isSigOdor & is_nontarget & upCountT>0,:);        
        bar(mean(dat1,1),'FaceColor','none','EdgeColor','b')%,'Normalization','probability')
        bar(mean(dat2,1),'FaceColor','none','EdgeColor','k')%,'Normalization','probability')                
        p2 = anovan([dat1(:);dat2(:)],[[ones(size(dat1(:)));2*ones(size(dat2(:)))] [repmat([1:size(dat1,2)]',[size(dat1,1) 1]);repmat([1:size(dat2,2)]',[size(dat2,1) 1])]],'display','off','model','interaction');
        ylabel('Proportion of cells')
        xlabel('Number of bins sig Target/NT')
        legend({'Target preference'; 'Non-target preference'},'Location','best')
        title(['All odor sig neurons, p = ' num2str(round(p2(3),2,'significant'))])

        subplot(5,2,10); hold on
        dat1 = upHistT(isUp & is_target & upCountT>0,:);
        dat2 = upHistT(isUp & is_nontarget & upCountT>0,:);
         bar(mean(dat1,1),'FaceColor','none','EdgeColor','b')%,'Normalization','probability')
        bar(mean(dat2,1),'FaceColor','none','EdgeColor','k')%,'Normalization','probability')        
        p2 = anovan([dat1(:);dat2(:)],[[ones(size(dat1(:)));2*ones(size(dat2(:)))] [repmat([1:size(dat1,2)]',[size(dat1,1) 1]);repmat([1:size(dat2,2)]',[size(dat2,1) 1])]],'display','off','model','interaction');
        ylabel('Proportion of cells')
        xlabel('Number of bins sig Target/NT')
        legend({'Target preference'; 'Non-target preference'},'Location','best')
        title(['Up, p = ' num2str(round(p2(3),2,'significant'))])
        
        set(gcf,'Position',[ 2019        -334         910        1308])
        helper_saveandclosefig([thisdir 'BriefLong_TargetMod_pvalue' num2str(sp.p_cutoff) '_' num2str(sum(times)*.2) '_Seconds_medcut' num2str(sp.medcut) '_' sp.ana_labs{iana}])
    end

    
    if sum(ismember(figs_to_run,0))>0
        % Chi-squared test between up and down modulated and modulated by target at all
        % Observed data
        n1 = sum(is_target & isOdorOff & isUp); N1 = sum(isUp & is_target);
        n2 = sum(is_nontarget & isOdorOff & isUp); N2 = sum(isUp & is_nontarget);
        % Pooled estimate of proportion
        p0 = (n1+n2) / (N1+N2);
        % Expected counts under H0 (null hypothesis)
        n10 = N1 * p0;
        n20 = N2 * p0;
        % Chi-square test, by hand
        observed = [n1 N1-n1 n2 N2-n2];
        expected = [n10 N1-n10 n20 N2-n20];
        [~,p0] = chi2gof([1 2 3 4],'freq',observed,'expected',expected,'ctrs',[1 2 3 4],'nparams',2);

        Utarget = [sum(is_target & ~isOdorOff & isUp) sum(is_target & isOdorOff & isUp)]./(sum(isUp & is_target));
        UNtarget = [sum(is_nontarget & ~isOdorOff & isUp) sum(is_nontarget & isOdorOff & isUp)]./(sum(isUp & is_nontarget));

        figure; hold on
        % subplot(2,2,[1 2]); hold on
        bar([Utarget; UNtarget],'stacked')
        % p01 = 1-binocdf(sum((is_target | is_nontarget) & isUp)-1,sum(isUp),sp.p_cutoff*2);
        % p02 = 1-binocdf(sum((is_target | is_nontarget) & isDown)-1,sum(isDown),sp.p_cutoff*2);
        set(gca,'xtick',1:2,'xticklabel',{'Target';'NonTarget'})
        text(-2.2,.2,{'Chi-squared test of Target vs NT'; 'being odor-off modulated';['p = ' num2str(round(p0,2,'significant'))]})
        % text(-2.2,.9,{'Binomial probability of Up';'cells being modulated,';['p = ' num2str(round(p01,2,'significant'))]})
        % text(-2.2,.6,{'Binomial probability of Down';'cells being modulated,';['p = ' num2str(round(p02,2,'significant'))]})
        xlim([-2.5 3])
        legend({'No Odor Off','OdorOff'},'Location','northeast')
        ylabel('Proportion of cells')

        helper_saveandclosefig([thisdir 'OdorOff_TargetMod_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana}])

        
        figure; hold on
        set(gcf,'Position',[2112         163         712         902])

        bb= subplot(4,2,1); hold on
        n = isOdorOff  & isUp; 
        dat = od(n,:);
        tbs = 1:size(dat,2);
        m1 = mean(zscore(dat,[],2)); sem1 = std(zscore(dat,[],2))./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
        rev = [m+sem]; fwd = [m-sem];
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
        plot(nanmean(zscore(dat,[],2)),'k');
        set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
        ylabel('Mean z-scored FR')
        title('Odor-Off and Odor-on Up')
        yl0 = get(gca,'ylim');

        aa = subplot(4,2,3); hold on
        m1 = mean(dat); sem1 = std(dat)./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
        rev = [m+sem]; fwd = [m-sem];
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
        plot(nanmean(dat),'k');
        set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
        ylabel('Mean FR')
        yl = get(gca,'ylim');
        title('Odor-Off and Odor-on Up')

        subplot(4,2,2); hold on
        n = ~isOdorOff & isUp; 
        dat = od(n,:);
        tbs = 1:size(dat,2);
        m1 = mean(zscore(dat,[],2)); sem1 = std(zscore(dat,[],2))./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
        rev = [m+sem]; fwd = [m-sem];
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
        plot(nanmean(zscore(dat,[],2)),'k');
        set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
        ylabel('Mean z-scored FR')
        y02 = get(gca,'ylim');
        set(gca,'ylim',[min([yl0(1) y02(1)]) max([yl0(2) y02(2)])])
        set(bb,'ylim',[min([yl0(1) y02(1)]) max([yl0(2) y02(2)])])
        title('Odor-on Up (not odor-off)')

        subplot(4,2,4); hold on
        m1 = mean(dat); sem1 = std(dat)./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
        rev = [m+sem]; fwd = [m-sem];
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
        plot(nanmean(dat),'k');
        set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
        ylabel('Mean FR')
        yl2 = get(gca,'ylim');
        set(gca,'ylim',[min([yl(1) yl2(1)]) max([yl(2) yl2(2)])])
        set(aa,'ylim',[min([yl(1) yl2(1)]) max([yl(2) yl2(2)])])
        title('Odor-on Up (not odor-off)')

        baselinebins = nanmean(squeeze(singledat.odorsel(:,sp.odorwindowlabel(:,2)<=0,1,3)),2);
        % figure; hold on; 
        subplot(4,2,[5:8]); hold on
        % histogram(baselinebins(~is_target & ~is_nontarget & isUp),0:.5:20,'Normalization','probability');
        histogram(baselinebins(isOdorOff & isUp),0:.5:20,'Normalization','probability'); %0:.2:1+max(baselinebins(isUp & (is_nontarget| is_target))))
        histogram(baselinebins(~isOdorOff  & isUp),0:.5:20,'Normalization','probability'); %0:.2:1+max(baselinebins(isUp & (is_nontarget| is_target))))
        legend('Odor-off','No odor-off')
        p = ranksum(baselinebins(~isOdorOff & isUp),baselinebins(isOdorOff & isUp),'tail','right')*2; % if sig, OdorOff<
        ylabel('Probability (of neuron)')
        xlabel('Average FR in 1 sec before odor')
        title(['Rank-sum that OdorOff is lower FR, p = ' num2str(round(p,2,'significant'))])
        helper_saveandclosefig([thisdir 'OdorOff_FRComp_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana}])




        figure; hold on
        set(gcf,'Position',[2112         163         712         902])

        subplot(4,2,1); hold on
        n = is_target & isUp; 
        dat = od(n,:);
        tbs = 1:size(dat,2);
        m1 = mean(zscore(dat,[],2)); sem1 = std(zscore(dat,[],2))./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
        rev = [m+sem]; fwd = [m-sem];
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
        plot(nanmean(zscore(dat,[],2)),'k');
        set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
        ylabel('Mean z-scored FR')
        title('Target preferring neurons')

        aa = subplot(4,2,3); hold on
        m1 = mean(dat); sem1 = std(dat)./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
        rev = [m+sem]; fwd = [m-sem];
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
        plot(nanmean(dat),'k');
        set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
        ylabel('Mean FR')
        yl = get(gca,'ylim');
        title('Target preferring neurons')

        subplot(4,2,2); hold on
        n = is_nontarget & isUp; 
        dat = od(n,:);
        tbs = 1:size(dat,2);
        m1 = mean(zscore(dat,[],2)); sem1 = std(zscore(dat,[],2))./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
        rev = [m+sem]; fwd = [m-sem];
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
        plot(nanmean(zscore(dat,[],2)),'k');
        set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
        ylabel('Mean z-scored FR')
        title('Non-target preferring neurons')

        subplot(4,2,4); hold on
        m1 = mean(dat); sem1 = std(dat)./sqrt(size(dat,1));
        m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
        rev = [m+sem]; fwd = [m-sem];
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
        plot(nanmean(dat),'k');
        set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
        ylabel('Mean FR')
        yl2 = get(gca,'ylim');
        set(gca,'ylim',[min([yl(1) yl2(1)]) max([yl(2) yl2(2)])])
        title('Non-target preferring neurons')

        set(aa,'ylim',[min([yl(1) yl2(1)]) max([yl(2) yl2(2)])])

        baselinebins = nanmean(squeeze(singledat.odorsel(:,sp.odorwindowlabel(:,2)<=0,1,3)),2);
        % figure; hold on; 
        subplot(4,2,[5:8]); hold on
        histogram(baselinebins(~is_target & ~is_nontarget & isUp),0:.5:20,'Normalization','probability');
        histogram(baselinebins(is_nontarget & isUp),0:.5:20,'Normalization','probability'); %0:.2:1+max(baselinebins(isUp & (is_nontarget| is_target))))
        histogram(baselinebins(is_target & isUp),0:.5:20,'Normalization','probability'); %0:.2:1+max(baselinebins(isUp & (is_nontarget| is_target))))
        legend('Not Modulated','Non-target','Target')
        p = ranksum(baselinebins(is_nontarget & isUp),baselinebins(is_target & isUp),'tail','left')*2; % if sig, Target>NonTarget
        ylabel('Probability (of neuron)')
        xlabel('Average FR in 1 sec before odor')
        title(['Rank-sum that Target>NonTarget, p = ' num2str(round(p,2,'significant'))])
        helper_saveandclosefig([thisdir 'TargetNT_FR_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana}])
    end

end
