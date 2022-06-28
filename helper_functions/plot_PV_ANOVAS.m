function plot_PV_ANOVAS(corr_repeat_trial_session,all_labels,dirs,itype)

thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'PV_ANOVAS/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

for itarget = 1:3
    if itarget<3
        dat = corr_repeat_trial_session(corr_repeat_trial_session(:,5)==itarget,:);
    else
        dat = corr_repeat_trial_session(corr_repeat_trial_session(:,5)>2,:);
    end

    for itt = 1:length(all_labels.ittlabs)
        [~,tbl] = anovan(dat(:,itt),dat(:,5:9),'varnames',{'RepeatNum','Trial','Session','Mouse','Rewarded'},'continuous',[2,3],'display','off');
        writecell(tbl,[thisdir all_labels.typelab{itype} '_' all_labels.labss{itarget} '_Anovan_table_' all_labels.ittlabs{itt} '_' all_labels.addon '.txt'],'Delimiter','tab')    
    end

    if itarget==3
        for itt = 2:3
            %self is bigger than other
            [p] = signrank(dat(:,1),dat(:,itt),'tail','right');

            %making sure once you factor in all these aspects that is still true
            dat1 = [dat(:,1);dat(:,itt)];     
            dat2 = [[ones(size(dat1,1)/2,1); ones(size(dat1,1)/2,1)*2] repmat(dat(:,5:9),[2 1])];
            [~,tbl,stats] = anovan(dat1,dat2,'varnames',{'Self/Other','RepeatNum','Trial','Session','Mouse','Rewarded'},'continuous',[3,4],'display','off');
            writecell(tbl,[thisdir all_labels.typelab{itype} '_PVPlot_Anovan_table_selfother' all_labels.addon '.txt'],'Delimiter','tab')    
            [~,m]  = multcompare(stats,'display','off');

            figure; hold on
            text(.2,.9,['Signrank p = ' num2str(p)])    
            text(.4,.8,['Self, median: ' num2str(median(dat(:,1),'omitnan'))])
            text(.4,.7,[all_labels.ittlabs{itt} ' , median: ' num2str(median(dat(:,itt),'omitnan'))])

            text(.2,.5,['Anovan, multiple comparison p = ' num2str(m(1,2))])    
            text(.4,.4,['Self, marginal mean: ' num2str(m(1,1))])
            text(.4,.3,[all_labels.ittlabs{itt} ' , marginal mean: ' num2str(m(2,1))]) 
            axis off
            helper_saveandclosefig([thisdir all_labels.typelab{itype} '_PVPlot_Anovan_SelfVOther_'  all_labels.ittlabs{itt} '_RewardedOnly_' all_labels.addon])
        end
    end        

end