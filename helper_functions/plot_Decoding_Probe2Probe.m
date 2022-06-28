function plot_Decoding_Probe2Probe(decode_fold_label,decode_fold_label_cat,all_labels,dirs,decoded_fold_group_all,itype)

thisdir = [dirs.figdir '/PopulationVectorsCorrelation/' ...
        'Probe2Probe/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

iplottype = 2;
icat = 1;
icor = 1;
for iraw = 1:2
    if iplottype == 1
      plotind = 4; plotlabel = 'Trial'; plotlabeloth = 'Session'; plotind2 = 5; binsz = 10;
    elseif iplottype == 2
      plotind = 5; plotlabel = 'Session'; plotlabeloth = 'Trial'; plotind2 = 4; binsz = 1;
    end

    if icat == 1
      OGdat = decode_fold_label; catlabel = 'Individual';
    elseif icat == 2
      OGdat = decode_fold_label_cat; catlabel = 'Categorical';
    end

    itt = 3;
    indall = decode_fold_label(:,1)>2 & decode_fold_label(:,2)==0;

     if icor ==1
         corlab = 'AllTrials';
     elseif icor == 2
         corlab = 'CorrectOnly';
         indall = indall & OGdat(:,7)==1;
     elseif icor == 3             
         corlab = 'IncorrectOnly';
         indall = indall & OGdat(:,7)==0;
     end

     dat = OGdat(indall,:);
     datS = decoded_fold_group_all(indall,3); % is it the other probes

    if iraw == 1
        datuse = datS; rawlab = 'Raw';
    elseif iraw==2
           [~,~,stats] = anovan(datS,dat(:,[plotind2 6:7]),'varnames',{plotlabeloth,'Mouse','Rewarded'},'continuous',[1],'display','off');    
            datuse = stats.resid; rawlab = 'Residuals';
    end
%         subplot(3,1,itt); hold on;
    %                     if iplottype==1

    figure; hold on
        plot(dat(:,plotind),datuse,'o','Color',[.5 .5 .5])
    %                     end
    dat1 = ceil((dat(:,plotind)-min(dat(:,plotind))+.001)/binsz);
    if sum(diff(sort(unique(dat1)))>1)==0
        m = splitapply(@mean,datuse,dat1);
        datss = splitapply(@std,datuse,dat1);
        h = hist(dat1,1:length(m));
        sem = datss./sqrt(h');
        rev = [m+sem]'; fwd = [m-sem]';
        tbs = [1:length(m)]*binsz;
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'blue');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 1];
        plot(tbs,m,'b');
    else
    %                         if iplottype==2
    %                              plot(dat(:,plotind),datuse,'o','Color',[.5 .5 .5])
    %                         end
    end
    [r,p] = corr(dat(:,plotind),datuse,'rows','complete');        
    tt = title([all_labels.ittlabs{itt+1} ', R = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))]);
    if p<.05 && r>0
        tt.Color = 'r';
    elseif p<.05 && r<0
        tt.Color = 'b';
    end
     xlabel(plotlabel)
     ylabel('Prob another probe')

     [~,tbl,~] = anovan(datS,dat(:,[4:7]),'varnames',{'Trial','Session','Mouse','Rewarded'},'continuous',[1],'display','off');  
    suptitle({[all_labels.typelab{itype} ' Decoding ' catlabel];['Over ' plotlabel 's, ' rawlab ', ' corlab];['Anovan, Session P = ' num2str(round(tbl{3,7},2,'significant'))]})
    set(gcf,'Position',[ 2128          80         533         829])
    helper_saveandclosefig([thisdir all_labels.typelab{itype} '_Decoding_' catlabel '_Over' plotlabel '_' rawlab '_' corlab all_labels.addon '_ProbeAsOtherProbe'])


    
     datS = decoded_fold_group_all(indall,2); % is it the NT

    if iraw == 1
        datuse = datS; rawlab = 'Raw';
    elseif iraw==2
           [~,~,stats] = anovan(datS,dat(:,[plotind2 6:7]),'varnames',{plotlabeloth,'Mouse','Rewarded'},'continuous',[1],'display','off');    
            datuse = stats.resid; rawlab = 'Residuals';
    end
%         subplot(3,1,itt); hold on;
    %                     if iplottype==1

    figure; hold on
        plot(dat(:,plotind),datuse,'o','Color',[.5 .5 .5])
    %                     end
    dat1 = ceil((dat(:,plotind)-min(dat(:,plotind))+.001)/binsz);
    if sum(diff(sort(unique(dat1)))>1)==0
        m = splitapply(@mean,datuse,dat1);
        datss = splitapply(@std,datuse,dat1);
        h = hist(dat1,1:length(m));
        sem = datss./sqrt(h');
        rev = [m+sem]'; fwd = [m-sem]';
        tbs = [1:length(m)]*binsz;
        p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'blue');
        p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 1];
        plot(tbs,m,'b');
    else
    %                         if iplottype==2
    %                              plot(dat(:,plotind),datuse,'o','Color',[.5 .5 .5])
    %                         end
    end

    [r,p] = corr(dat(:,plotind),datuse,'rows','complete');        
    tt = title([all_labels.ittlabs{itt+1} ', R = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))]);
    if p<.05 && r>0
        tt.Color = 'r';
    elseif p<.05 && r<0
        tt.Color = 'b';
    end
     xlabel(plotlabel)
     ylabel('Prob NT')
     [~,tbl,~] = anovan(datS,dat(:,[4:7]),'varnames',{'Trial','Session','Mouse','Rewarded'},'continuous',[1],'display','off');  
    suptitle({[all_labels.typelab{itype} ' Decoding ' catlabel];['Over ' plotlabel 's, ' rawlab ', ' corlab];['Anovan, Session P = ' num2str(round(tbl{3,7},2,'significant'))]})
    set(gcf,'Position',[ 2128          80         533         829])
    helper_saveandclosefig([thisdir all_labels.typelab{itype} '_Decoding_' catlabel '_Over' plotlabel '_' rawlab '_' corlab all_labels.addon '_ProbeAsNT'])
end