function plot_Selectivity_OdorOff(singledat,selectivity_dat,selectivity_params,dirs,figs_to_run)
% one panel for Figure s5

od = singledat.toplotOdorLick(:,1:6000);


sp = selectivity_params;
sd = selectivity_dat;

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'OdorOff/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

  %%%%% Odor-off responses
lablabs = {'Odor on positive neurons','Odor on negative neurons', 'Odor on positive neurons and odor off',...
    'Odor on negative neurons and odor off','Odor-off only'};

if sum(ismember(figs_to_run,[0 15]))>0
    for iana = 1        
        alldat = cat(2,sd.all_Neurons(:,3:4), sd.all_Neurons(:,3) & sd.all_Neurons(:,13), sd.all_Neurons(:,4) & sd.all_Neurons(:,13),~sd.all_Neurons(:,3) & ~sd.all_Neurons(:,4) & sd.all_Neurons(:,13));        
        
        
         ntotal = size(sd.neuron_info,1); 
        
        if iana==2    
            ntotal = sum(sd.neuron_info(:,7));
            alldat = alldat & repmat(sd.neuron_info(:,7)==1,[1 size(alldat,2)]);
        end
        
        figure; hold on
        for idat = 1:5
            subplot(5,2,(idat-1)*2+1); hold on
            n = find(alldat(:,idat));             
            title([lablabs{idat} ', ' num2str(length(n))]) %round(100*length(n)./length(alldat(:,idat)),2,'significant')) '%'])
            if length(n)>1
                dat = od(n,:);

                dat2 = zscore(dat,[],2);
                datBS = dat2-mean(dat2(:,1:1000),2);

                [~,mm] = max(dat(:,3000:-1:1000),[],2);
                [~,ord] = sort(mm);
                imagesc(dat(ord,:)./max(dat(ord,:),[],2))
                set(gca,'ytick',1:size(dat,1),'yticklabel',n(ord))
                set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
                colorbar
                axis tight ij
                
                subplot(5,2,(idat-1)*2+2); hold on
                tbs = 1:size(dat,2);
                m1 = mean(datBS); sem1 = std(datBS)./sqrt(size(dat,1));
                m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
                rev = [m+sem]; fwd = [m-sem];
                p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
                p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
                plot(nanmean(datBS),'k');
                set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
            end
        end
        
        set(gcf,'Position',[1958        -307        1003        1633])
        helper_saveandclosefig([thisdir 'Neurons' num2str(ntotal) '_pvalue' num2str(sp.p_cutoff) '_' sp.ana_labs{iana}])
    end      
end


if sum(ismember(figs_to_run,0))>0
    baselinebins = nanmean(squeeze(singledat.odorsel(:,sp.odorwindowlabel(:,2)<=2000,1,4)),2);
    lablabs = {'Odor-off & up only';'Odor-off,up only and T/NT mod';};
    for iana = 1        
    
        alldat = cat(2,sd.all_Neurons(:,14) & ~sd.all_Neurons(:,3) & ~sd.all_Neurons(:,4), ...
            sd.all_Neurons(:,14) & ~sd.all_Neurons(:,3) & ~sd.all_Neurons(:,4) & (sd.all_Neurons(:,5) | ~sd.all_Neurons(:,6)),...
            sd.all_Neurons(:,15), sd.all_Neurons(:,16));        
        
         ntotal = size(sd.neuron_info,1); 
        
        if iana==2    
            ntotal = sum(sd.neuron_info(:,7));
            alldat = alldat & repmat(sd.neuron_info(:,7)==1,[1 size(alldat,2)]);
        end
        
        figure; hold on;
        for idat = 1:2            
            subplot(3,2,(idat-1)*2+1); hold on
            n = find(alldat(:,idat)); 
            dat = od(n,:);

            dat2 = zscore(dat,[],2);
            datBS = dat2-mean(dat2(:,1:1000),2);

            [~,mm] = max(dat(:,4000:-1:1000),[],2);
            [~,ord] = sort(mm);
            imagesc(dat(ord,:)./max(dat(ord,:),[],2))
            set(gca,'ytick',1:size(dat,1),'yticklabel',n(ord))
            set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
            colorbar
            axis tight ij
            title([ lablabs{idat} ', ' num2str(length(n))]) % round(100*length(n)./ntotal,2,'significant')) '%'])
            xlabel('Time since odor-on')
    
            if size(dat,1)>1
                subplot(3,2,(idat-1)*2+2); hold on
                tbs = 1:size(dat,2);
                m1 = mean(datBS); sem1 = std(datBS)./sqrt(size(dat,1));
                m = smoothdata(m1,'movmean',sp.sm); sem = smoothdata(sem1,'movmean',sp.sm);
                rev = [m+sem]; fwd = [m-sem];
                p1 = patch([tbs tbs(length(tbs):-1:1)],[fwd rev(end:-1:1)],'black');
                p1.FaceAlpha=.2; p1.EdgeAlpha=0;    p1.FaceColor = [0 0 0];
                plot(nanmean(datBS),'k');
                set(gca,'xtick',[1000:1000:6000],'xticklabel',[0:5])
                xlabel('Time since odor-on')
            end
        end
    
    
        n = find(alldat(:,1));
       
        if length(n)<2
            continue
        end
        
        subplot(3,2,5); hold on
        histogram(baselinebins(n),0:.2:max(baselinebins(n)))
        title(['Mean: ' num2str(round(mean(baselinebins(n)),2,'significant')) ',Median: ' num2str(round(median(baselinebins(n)),2,'significant'))])
        xlabel('Baseline (-1-2sec) FR (.2 bins)')
        ylabel('Number of neurons')
    
        subplot(3,2,6); hold on
        [~,ord] = sort(baselinebins(n));
        % yyaxis left
        barh(1:length(n),baselinebins(n(ord)))
        set(gca,'ytick',1:length(n),'yticklabel',n(ord))
        n2 = ismember(n,find(alldat(:,3))); 
        barh(find(n2),baselinebins(n(ord(n2))))
        n2 = ismember(n,find(alldat(:,4))); 
        barh(find(n2),baselinebins(n(ord(n2))))
        legend('Not T/NT mod','T mod','NT mod')
        xlabel('Baseline (-1-2sec) FR')
    
        set(gcf,'Position',[1921        -350        1080        1803])
    
        helper_saveandclosefig([thisdir 'Neurons' num2str(ntotal) '_pvalue' num2str(sp.p_cutoff) '_VenkiNeurons_FR_' sp.ana_labs{iana}])
    end
end

close all