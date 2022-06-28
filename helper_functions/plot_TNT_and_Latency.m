function plot_TNT_and_Latency(corr_repeat_trial_session,all_labels,dirs,itype,across_lookback,figs_to_run)

numshuff = 1000;
lick_cutoff_ms = 500;

thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'TNT_and_Latency/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

% corr_repeat_trial_session
% 1. PVcorr self, 2. PVcorr target, 3. PVcorr nontarget, 4. PVcorrOthers, 
% 5. trialtype#, 6. trial, 7. session, 8. mouse, 9. accuracy, 
% 10. lookback, 11. lick latency

dat2 = corr_repeat_trial_session(corr_repeat_trial_session(:,5)>2,:);
delind = dat2(:,11)<lick_cutoff_ms;
dat2(delind,:) = [];
pbs = unique(dat2(:,5));

%testing if true for T/NT, it is
% datNT = corr_repeat_trial_session(corr_repeat_trial_session(:,5)<3,:);
% delind2 = datNT(:,11)<lick_cutoff_ms;
% datNT(delind2,:) = [];

ALB = across_lookback(corr_repeat_trial_session(:,5)>2,:);
ALB(delind,:) = [];
% all_labels.LB
% plot across lookback paramters
Ps = NaN(length(all_labels.LB),1);
for iLB = 1:length(all_labels.LB)
    Ps(iLB,1) = ranksum(ALB(dat2(:,9)==1,iLB),ALB(dat2(:,9)==0,iLB),'tail','left');
end

figure; hold on
set(gcf,'Position',[ 2362         284         455         388])
plot(all_labels.LB,Ps,'ok')
xl = [min(all_labels.LB)-1 max(all_labels.LB)+1];
plot(xl,[0.05 0.05],'--r')
set(gca,'xlim',xl)
ylabel('P-value')
xlabel('Number of trials back to average over')
helper_saveandclosefig([thisdir 'LatencyVsPVcorr_Lookback_' all_labels.typelab{itype} all_labels.addon])

%Correct trials tend to have longer latencies
for iprobe = 1:length(pbs)+1
    
    if iprobe<=length(pbs) && sum(ismember(figs_to_run,0))==0
        continue
    end

    if iprobe > length(pbs)
        dat = dat2;
    else
        dat = dat2(dat2(:,5)==pbs(iprobe),:);
    end

    figure; hold on; 
    set(gcf,'Position',[ 680    90   997   888])
    
    subplot(2,2,1); hold on; 
    plot(dat(dat(:,9)==1,11),dat(dat(:,9)==1,10),'ob','MarkerSize',2); plot(dat(dat(:,9)==0,11),dat(dat(:,9)==0,10),'xr')
    xl = get(gca,'xlim'); yl = get(gca,'ylim');
    [r1,p1] = corr(dat(dat(:,9)==1,11),dat(dat(:,9)==1,10),'rows','complete');
    [r0,p0] = corr(dat(dat(:,9)==0,11),dat(dat(:,9)==0,10),'rows','complete');
    [r,p] = corr(dat(:,11),dat(:,10),'rows','complete');    
    rps = NaN(numshuff,2);
    for ishuff = 1:numshuff
        [rps(ishuff,1),rps(ishuff,2) ] = corr(dat(randperm(size(dat,1)),11),dat(randperm(size(dat,1)),10),'rows','complete');
    end
    Pshuff = (1+sum(abs(rps(:,1))>abs(r)))./(1+numshuff); %two-sided p = 9.99e-4
    Pshuff = (1+sum(rps(:,1)<r))./(1+numshuff); %one-sided p = 9.99e-4

    text(xl(1)+range(xl)/4,yl(1)+(range(yl)/40)*10,['r = ' num2str(round(r0,2,'significant')) ', p = ' num2str(round(p0,2,'significant'))],'Color','r')    
    text(xl(1)+range(xl)/4,yl(1)+(range(yl)/40)*8,['r = ' num2str(round(r1,2,'significant')) ', p = ' num2str(round(p1,2,'significant'))],'Color','b')
    text(xl(1)+range(xl)/4,yl(1)+(range(yl)/40)*6,['n = ' num2str(sum(~isnan(dat(:,11)) & ~isnan(dat(:,10)))) ', r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))])
    text(xl(1)+range(xl)/4,yl(1)+(range(yl)/40)*4,['Shuffp' num2str(numshuff) ' = ' num2str(round(Pshuff,2,'significant'))])
    ylabel('PV correlation T-NT')
    xlabel('Latency (ms)')    
    legend('Correct','Incorrect')
    
    subplot(2,2,3); hold on
    histogram(dat(dat(:,9)==1,11),xl(1):50:xl(2))
    histogram(dat(dat(:,9)==0,11),xl(1):50:xl(2))
    xlabel('Latency (ms)')
    ylabel('Count')    
    yl2 = get(gca,'ylim');
    plot([median(dat(dat(:,9)==1,11),'omitnan') median(dat(dat(:,9)==1,11),'omitnan')],yl2,'b-')
    plot([median(dat(dat(:,9)==0,11),'omitnan') median(dat(dat(:,9)==0,11),'omitnan')],yl2,'r-')
    I1 = dat(dat(:,9)==1,11);
    I2 = dat(dat(:,9)==0,11);
    P = ranksum(I1,I2,'tail','right');
    text(1100,yl2(2)-range(yl2)/2,['Nc = ' num2str(sum(~isnan(I1))) ', Nin = ' num2str(sum(~isnan(I2))) ', P = ' num2str(round(P,3,'significant'))])
    legend('Correct','Incorrect','Correct','Incorrect','Location','Best')
    
    
    subplot(2,2,2); hold on; 
    histogram(dat(dat(:,9)==1,10),(yl(1):.01:yl(2)),'Orientation','horizontal'); %,'Normalization','probability');
    histogram(dat(dat(:,9)==0,10),(yl(1):.01:yl(2)),'Orientation','horizontal'); %,'Normalization','probability');
    ylabel('PV correlation T-NT')
    xlabel('Count')
    xl2 = get(gca,'xlim'); yl2 = get(gca,'ylim');
    plot(xl2,[median(dat(dat(:,9)==1,10),'omitnan') median(dat(dat(:,9)==1,10),'omitnan')],'b-')
    plot(xl2,[median(dat(dat(:,9)==0,10),'omitnan') median(dat(dat(:,9)==0,10),'omitnan')],'r-')
    legend('Correct','Incorrect','Correct','Incorrect')
    I1 = dat(dat(:,9)==1,10);
    I2 = dat(dat(:,9)==0,10);
    P = ranksum(I1,I2,'tail','left');
    text(xl2(1)+range(xl2)/4,yl2(1)+range(yl2)/4,['Nc = ' num2str(sum(~isnan(I1))) ', Nin = ' num2str(sum(~isnan(I2))) ', P = ' num2str(round(P,3,'significant'))])

    %correct-incorrect
    r = median(I1)-median(I2);
    rps = NaN(numshuff,1);
    for ishuff = 1:numshuff
        shufcor = dat(randperm(size(dat,1)),9);
        rps(ishuff,1) = median(dat(shufcor==1,10))...
            -median(dat(shufcor==0,10));        
    end
    Pshuff = (1+sum(rps(:,1)<r))./(1+numshuff); %one-sided, p = 0.0155
%     Pshuff = (1+sum(abs(rps(:,1))>abs(r)))./(1+numshuff); %two-sided, p =0.0318
 
    text(xl2(1)+range(xl2)/4,yl2(1)+range(yl2)/3,['Shuffp' num2str(numshuff) ' = ' num2str(round(Pshuff,3,'significant'))])
 
    if iprobe > length(pbs)
        helper_saveandclosefig([thisdir 'LatencyVsPVcorr_' all_labels.typelab{itype} all_labels.addon])
    else
        helper_saveandclosefig([thisdir 'LatencyVsPVcorr_' all_labels.typelab{itype} all_labels.addon '_probe ' num2str(iprobe)])
    end
end




if sum(ismember(figs_to_run,0))>0
    dat = dat2;
    
      %What correlates with the T/NT distance 
    % (latency does even with other factors, accuracy does not, not as strong)
    [~,tbl] = anovan(dat(:,10),dat(:,[6:9 11]),'varnames',...
        {'Trial','Session','Mouse','Accuracy','Latency'},...
        'continuous',[1:2 5],'display','off');    
    writecell(tbl,[thisdir 'PredictTNT_' ...
        all_labels.typelab{itype} '_ANOVAN.txt'],'Delimiter','tab')    
    
    [~,tbl] = anovan(dat(:,10),dat(:,[6:9 11]),'varnames',...
        {'Trial','Session','Mouse','Accuracy','Latency'},...
        'continuous',[1:2 5],'model','interaction','display','off');    
    writecell(tbl,[thisdir 'PredictTNT_' ...
        all_labels.typelab{itype} '_ANOVAN_interaction.txt'],...
        'Delimiter','tab')        
    

    %What correlates with the Accuracy (latency does, not TNT) 
    [~,tbl] = anovan(dat(:,9),dat(:,[6:8 10:11]),'varnames',...
        {'Trial','Session','Mouse','TNT','Latency'},...
        'continuous',[1:2 4:5],'display','off');    
    writecell(tbl,[thisdir 'PredictAcc_' ...
        all_labels.typelab{itype} '_ANOVAN.txt'],'Delimiter','tab')    

    % TNT does here - only model with significant TNT/accuracy
    [~,tbl] = anovan(dat(:,9),dat(:,[6:8 10:11]),'varnames',...
        {'Trial','Session','Mouse','TNT','Latency'},...
        'continuous',[1:2 4:5],'model','interaction','display','off');
    writecell(tbl,[thisdir 'PredictAcc_' ...
        all_labels.typelab{itype} '_ANOVAN_interaction.txt'],...
        'Delimiter','tab')    


    %What correlates with the T/NT distance 
    % (latency does even with other factors, accuracy does not, not as strong)
    [~,tbl] = anovan(dat(:,10),dat(:,[5:9 11]),'varnames',...
        {'Probe','Trial','Session','Mouse','Accuracy','Latency'},...
        'continuous',[2:3 6],'display','off');    
    writecell(tbl,[thisdir 'PredictTNT_' ...
        all_labels.typelab{itype} '_ANOVAN_Probe.txt'],'Delimiter','tab')    
    
    [~,tbl] = anovan(dat(:,10),dat(:,[5:9 11]),'varnames',...
        {'Probe','Trial','Session','Mouse','Accuracy','Latency'},...
        'continuous',[2:3 6],'model','interaction','display','off');    
    writecell(tbl,[thisdir 'PredictTNT_' ...
        all_labels.typelab{itype} '_ANOVAN_Probe_interaction.txt'],...
        'Delimiter','tab')        
    

    %What correlates with the Accuracy (latency does, not TNT) 
    [~,tbl] = anovan(dat(:,9),dat(:,[5:8 10:11]),'varnames',...
        {'Probe','Trial','Session','Mouse','TNT','Latency'},...
        'continuous',[2:3 5:6],'display','off');    
    writecell(tbl,[thisdir 'PredictAcc_' ...
        all_labels.typelab{itype} '_ANOVAN_Probe.txt'],'Delimiter','tab')    

    [~,tbl] = anovan(dat(:,9),dat(:,[5:8 10:11]),'varnames',...
        {'Probe','Trial','Session','Mouse','TNT','Latency'},...
        'continuous',[2:3 5:6],'model','interaction','display','off');
    writecell(tbl,[thisdir 'PredictAcc_' ...
        all_labels.typelab{itype} '_ANOVAN_Probe_interaction.txt'],...
        'Delimiter','tab')    

    %What correlates with the lick latency
     [~,tbl] = anovan(dat(:,11),dat(:,[5:9 10]),'varnames',...
        {'Probe','Trial','Session','Mouse','Accuracy','TNT'},...
        'continuous',[2:3 6],'display','off');    
    writecell(tbl,[thisdir 'PredictLatency_' ...
        all_labels.typelab{itype} '_Probe_ANOVAN.txt'],'Delimiter','tab')    
    
    [~,tbl] = anovan(dat(:,11),dat(:,[5:9 10]),'varnames',...
        {'Probe','Trial','Session','Mouse','Accuracy','TNT'},...
        'continuous',[2:3 6],'model','interaction','display','off');    
    writecell(tbl,[thisdir 'PredictLatency_' ...
        all_labels.typelab{itype} '_Probe_ANOVAN_interaction.txt'],...
        'Delimiter','tab')        
    
    mouseses = NaN(max(dat(:,7)),max(dat(:,8)),2);
    figure; hold on
    for imouse = 1:max(dat(:,8))
        for ises = 1:max(dat(:,7))
            if sum(dat(:,7)==ises & dat(:,8)==imouse)>0
                ind = dat(:,7)==ises & dat(:,8)==imouse;
                mouseses(ises,imouse,1) = mean(dat(ind,9),'omitnan'); %accuracy on probe trials                
                mouseses(ises,imouse,2) = mean(dat(ind,10),'omitnan'); %T-NT correlation
            end
        end
        plot(mouseses(:,imouse,1),mouseses(:,imouse,2),'o','MarkerSize',20)
    end
    t1 = mouseses(:,:,1); t2 = mouseses(:,:,2); [r,p] = corr(t1(:),t2(:),'rows','complete');
    ylabel('PV correlation T-NT')
    xlabel('Accuracy on Probe trials')
    title(['r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))])
    helper_saveandclosefig([thisdir 'LatencyVsPVcorr_Session' all_labels.typelab{itype} all_labels.addon])
    
    
    % we are still getting better overall (during probes, not during NTR)
    dat = corr_repeat_trial_session;
    [~,dat1,stats1] = anovan(dat(:,9),dat(:,6:8),'varnames',{'Trial','Session','Mouse'},'continuous',1:2,'nested',[ 0 0 0; 0 0 1 ;0 0 0],'display','off');    
    
    % getting better at probes (not NTRs)
    dat = corr_repeat_trial_session(corr_repeat_trial_session(:,5)>2,:);
    [~,dat2,stats2] = anovan(dat(:,9),dat(:,6:8),'varnames',{'Trial','Session','Mouse'},'continuous',1:2,'nested',[ 0 0 0; 0 0 1 ;0 0 0],'display','off');    
    
    % getting better at targets
    dat = corr_repeat_trial_session(corr_repeat_trial_session(:,5)==1,:);
    [~,dat3,stats3] = anovan(dat(:,9),dat(:,6:8),'varnames',{'Trial','Session','Mouse'},'continuous',1:2,'nested',[ 0 0 0; 0 0 1 ;0 0 0],'display','off');
    
    % not getting better at NTs during Probes (are getting worse during NTRs)
    dat = corr_repeat_trial_session(corr_repeat_trial_session(:,5)==2,:);
    [~,dat4,stats4] = anovan(dat(:,9),dat(:,6:8),'varnames',{'Trial','Session','Mouse'},'continuous',1:2,'nested',[ 0 0 0; 0 0 1 ;0 0 0],'display','off');    
    
    figure; hold on;
    text(.01,.9,['Overall change in accuracy over sessions (' all_labels.typelab{itype} ' sessions):'])
    text(.03,.8,['Model: ' helper_text_rp(stats1.coeffs(3),dat1{3,7},2)]) % ', corr: ' helper_text_rp(r1,p1,2)])
    text(.01,.7,['Change in accuracy of ' all_labels.typelab{itype} ' over sessions:'])
    text(.03,.6,['Model: ' helper_text_rp(stats2.coeffs(3),dat2{3,7},2)]) % ', corr: ' helper_text_rp(r2,p2,2)])
    text(.01,.5,'Overall change in accuracy of target over sessions:')
    text(.03,.4,['Model: ' helper_text_rp(stats3.coeffs(3),dat3{3,7},2)]) % ', corr: ' helper_text_rp(r3,p3,2)])
    text(.01,.3,'Overall change in accuracy of nontarget over sessions: ')
    text(.03,.2,['Model: ' helper_text_rp(stats4.coeffs(3),dat4{3,7},2)]) % ', corr: ' helper_text_rp(r4,p4,2)])
    axis off
    helper_saveandclosefig([thisdir 'Acc_overSession_' all_labels.typelab{itype} all_labels.addon])
    
    
    
    %%%% Lick latencies on correct trials
    dat = corr_repeat_trial_session(corr_repeat_trial_session(:,9)==1,:);
    [~,dat1,stats1] = anovan(dat(:,11),dat(:,6:8),'varnames',{'Trial','Session','Mouse'},'continuous',1:2,'nested',[ 0 0 0; 0 0 1 ;0 0 0],'display','off');    
    
    % getting better at probes (not NTRs)
    dat = corr_repeat_trial_session(corr_repeat_trial_session(:,5)>2 & corr_repeat_trial_session(:,9)==1,:);
    [~,dat2,stats2] = anovan(dat(:,11),dat(:,6:8),'varnames',{'Trial','Session','Mouse'},'continuous',1:2,'nested',[ 0 0 0; 0 0 1 ;0 0 0],'display','off');
    
    dat = corr_repeat_trial_session(corr_repeat_trial_session(:,5)==1 & corr_repeat_trial_session(:,9)==1,:);
    [~,dat3,stats3] = anovan(dat(:,11),dat(:,6:8),'varnames',{'Trial','Session','Mouse'},'continuous',1:2,'nested',[ 0 0 0; 0 0 1 ;0 0 0],'display','off');
    
    dat = corr_repeat_trial_session(corr_repeat_trial_session(:,5)==2 & corr_repeat_trial_session(:,9)==1,:);
    [~,dat4,stats4] = anovan(dat(:,11),dat(:,6:8),'varnames',{'Trial','Session','Mouse'},'continuous',1:2,'nested',[ 0 0 0; 0 0 1 ;0 0 0],'display','off');    
    
    
    figure; hold on;
    text(.01,.9,['Overall change in lick latency (' all_labels.typelab{itype} ' sessions):'])
    text(.03,.8,['Model: ' helper_text_rp(stats1.coeffs(3),dat1{3,7},2)]) % ', corr: ' helper_text_rp(r1,p1,2)])
    text(.01,.7,['Change in lick latency of ' all_labels.typelab{itype} ' over sessions:'])
    text(.03,.6,['Model: ' helper_text_rp(stats2.coeffs(3),dat2{3,7},2)]) % ', corr: ' helper_text_rp(r2,p2,2)])
    text(.01,.5,'Overall change in lick latency of target over sessions:')
    text(.03,.4,['Model: ' helper_text_rp(stats3.coeffs(3),dat3{3,7},2)]) % ', corr: ' helper_text_rp(r3,p3,2)])
    text(.01,.3,'Overall change in lick latency of nontarget over sessions: ')
    text(.03,.2,['Model: ' helper_text_rp(stats4.coeffs(3),dat4{3,7},2)]) % ', corr: ' helper_text_rp(r4,p4,2)])
    axis off
    helper_saveandclosefig([thisdir 'Lick_overSession_' all_labels.typelab{itype} all_labels.addon])
end