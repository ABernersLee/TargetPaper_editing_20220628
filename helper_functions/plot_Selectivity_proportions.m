function plot_Selectivity_proportions(singledat,selectivity_dat,selectivity_params,dirs,figs_to_run)

%%% this needs to be cleaned up
% using 2 figs for figure 5

thisdir = [dirs.figdir '/TypesOfResponses/' ...
        'Selectivity_proportions/'];

if ~isfolder(thisdir)
    mkdir(thisdir)
end

sp = selectivity_params;
sd = selectivity_dat;

isProbeRepeatSession(:,:,1) = any(squeeze(sum(sum(~isnan(singledat.odorselective(:,:,:,:,1,15:17)),4),3))>0,3);
isProbeRepeatSession(:,:,2) = any(squeeze(sum(sum(~isnan(singledat.odorselective(:,:,:,:,1,18:19)),4),3))>0,3);


for itype = 1:2
    A1 = squeeze(singledat.odorselective(:,:,:,sp.OdorIndSmall,2,3));        
    S = sum(squeeze(any(any(~isnan(A1),3),4)),3);
    pcutoff = (sp.p_cutoff./(sum(sp.OdorIndSmall)*S));
    pcutoff(S==0) = NaN;        
    A11 = squeeze(sum(sum(any(A1<repmat(pcutoff,[1 1 size(A1,3)...
        size(A1,4) size(A1,5)]),4),3),5)); 
    indd = squeeze(sum(sum(any(~isnan(A1),4),3),5)); 
    sess = repmat([1:size(A11,2)],[size(A11,1) 1 size(A11,3)]);        
    ani = repmat([1:size(A11,1)],[size(A11,2) 1 size(A11,3)])';
    datA = A11./indd;   
    datA(isProbeRepeatSession(:,:,itype)==0) = NaN;
    [p2,~,~] = anovan(datA(:),[sess(:) ani(:)],'varnames',...
        {'Session';'Mouse'},'continuous',1,'display','off');
    
    if sum(ismember(figs_to_run,0))>0
        figure; hold on
        plot(sess(:),datA(:),'ko')
        [r,p] = corr(sess(:),datA(:),'rows','complete');
        N = sum(~isnan(sess(:)) & ~isnan(datA(:))); 
        suptitle({['Prop Corr Sig; N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{itype}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})         
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' sp.rep_probe{itype} ...
            '_AnyBins'])  
    end
    
    A1 = squeeze(singledat.odorselective_all(:,:,:,1,2));       
    A11 = sum(A1<sp.p_cutoff,3);
    indd = sum(~isnan(A1),3); %squeeze(sum(sum(any(~isnan(A1),4),3),5)); 
    sess = repmat([1:size(A11,2)],[size(A11,1) 1 size(A11,3)]);        
    ani = repmat([1:size(A11,1)],[size(A11,2) 1 size(A11,3)])';
    datA = A11./indd;    
    datA(isProbeRepeatSession(:,:,itype)==0) = NaN;    
    [p2,tbl2,~] = anovan(datA(:),[sess(:) ani(:)],'varnames', ...
        {'Session';'Mouse'},'continuous',1,'display','off');
        
    if sum(ismember(figs_to_run,0))>0
        figure; hold on
        plot(sess(:),datA(:),'ko')
        [r,p] = corr(sess(:),datA(:),'rows','complete');
        N = sum(~isnan(sess(:)) & ~isnan(datA(:))); 
        xlabel('Sessions'); ylabel('Proportion of neurons')
        suptitle({['Prop Corr Sig; N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{itype}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})   
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' sp.rep_probe{itype} ...
            '_OverallP'])    
    
        figure; hold on
        for imouse = 1:4
            if sum(~isnan(datA(imouse,:)))>0
                plot(sess(imouse,:),datA(imouse,:),'o');          
            end
        end    
        [r,p] = corr(sess(:),datA(:),'rows','complete');
        N = sum(~isnan(sess(:)) & ~isnan(datA(:))); 
        xlabel('Sessions'); ylabel('Proportion of neurons')
           suptitle({['Prop Corr Sig; N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{itype}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})         
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' sp.rep_probe{itype} ...
            '_OverallP_mouseplot'])    
    end

    %figure 5e
    if (sum(ismember(figs_to_run,[0 5]))>0 && itype==2) || (sum(ismember(figs_to_run,0))>0)
        figure; hold on
        for imouse = 1:4
            if sum(~isnan(datA(imouse,:)))>0
                p = plot(sess(imouse,:),datA(imouse,:),'o-');
    %             [xF,yF] = helper_plotbestfitline(sess(imouse,:),datA(imouse,:));
    %             plot(xF,yF,'Color',p.Color,'LineWidth',1)
            end
        end
        [xF,yF] = helper_plotbestfitline(sess(:),datA(:));
        plot(xF,yF,'Color','k','LineWidth',2)
        [r,p] = corr(sess(:),datA(:),'rows','complete');
        N = sum(~isnan(sess(:)) & ~isnan(datA(:))); 
        axis tight
        xl = get(gca,'xlim');
        set(gca,'xlim',[xl(1)-1 xl(2)+2])
        xlabel('Sessions'); ylabel('Proportion of Target/NT mod. neurons')
           suptitle({['Prop Corr Sig; N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{itype}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})         
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' sp.rep_probe{itype} ...
            '_OverallP_mouseplotline'])   
         writecell(tbl2,[thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' sp.rep_probe{itype} ...
            '_OverallP_mouseplotline_ANOVAN.txt'],'Delimiter','tab')  
         writecell(tbl2,[thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' sp.rep_probe{itype} ...
            '_OverallP_mouseplotline_ANOVAN.xlsx'])  
    end


    if sum(ismember(figs_to_run,0))>0
        sess2 = NaN(size(sess));
        for imouse = 1:4
            sess2(imouse,~isnan(datA(imouse,:))) = 1:sum(~isnan(datA(imouse,:)));
        end
    
         figure; hold on
        for imouse = 1:4
            if sum(~isnan(datA(imouse,:)))>0
                p = plot(sess2(imouse,:),datA(imouse,:),'o');
                [xF,yF] = helper_plotbestfitline(sess2(imouse,:),datA(imouse,:));
                plot(xF,yF,'Color',p.Color,'LineWidth',1)
            end
        end
        [xF,yF] = helper_plotbestfitline(sess2(:),datA(:));
        plot(xF,yF,'Color','k','LineWidth',2)
        [r,p] = corr(sess2(:),datA(:),'rows','complete');
        N = sum(~isnan(sess2(:)) & ~isnan(datA(:))); 
        xlabel('Sessions'); ylabel('Proportion of neurons')
           suptitle({['Prop Corr Sig; N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{itype}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})         
           axis tight
           yl = get(gca,'ylim');
           xlim([0 19])       
    %      ylim([yl(1)-.1 yl(2)+.1])
        ylim([0 .7])
         set(gca,'xtick',5:5:15,'ytick',0:.2:1)
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' sp.rep_probe{itype} ...
            '_OverallP_mouseplotline_SessionsInOrder'])    
    
        sess2 = NaN(size(sess));
        sess3 = NaN(size(sess));
        for imouse = 1:4
            sess2(imouse,1:sum(~isnan(datA(imouse,:)))) = datA(imouse,~isnan(datA(imouse,:)));
            sess3(imouse,1:sum(~isnan(datA(imouse,:)))) = 1:sum(~isnan(datA(imouse,:)));
        end
        figure; hold on
        for imouse = 1:4
            if sum(~isnan(sess2(imouse,:)))>0
                plot(sess2(imouse,:),'o-');
            end
        end
    %     errorbar(mean(sess2,1,'omitnan'),std(sess2,1,'omitnan')./sqrt(sum(~isnan(sess2),1)),'k','LineWidth',2)
       [xF,yF] = helper_plotbestfitline(1:sum(~isnan(mean(sess2,1,'omitnan'))),mean(sess2(:,1:sum(~isnan(mean(sess2,1,'omitnan')))),1,'omitnan'));
        plot(xF,yF,'Color','k','LineWidth',2)
        [r,p] = corr(sess3(:),sess2(:),'rows','complete');
        N = sum(~isnan(sess2(:)) & ~isnan(datA(:))); 
        xlabel('Sessions'); ylabel('Proportion of neurons')
           suptitle({['Prop Corr Sig; N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{itype}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})         
           axis tight
           yl = get(gca,'ylim');
           xlim([0 19])       
    %      ylim([yl(1)-.1 yl(2)+.1])
        ylim([0 .7])
         set(gca,'xtick',5:5:15,'ytick',0:.2:1)
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' ...
            sp.ana_labs{3} '_TargetNonTargt_PropCorrSig_' sp.rep_probe{itype} ...
            '_OverallP_mouseplotline_SessionsInOrder2'])    
    end

end



%%%%% paired comparison of trial downsampled target vs probe or repeat    
repnum = [5 7 10 12 15 17;8 9 13 14 18 19];
allp = NaN(size(singledat.odorselective,1),size(singledat.odorselective,2),2,2,3); 
%probe or repeat / target, whether probe or repeat, # of odor
for irep = 1:2

    allp2 = NaN(size(singledat.odorselective,1),size(singledat.odorselective,2),3,3); 

    A1 = squeeze(singledat.odorselective(:,:,:,sp.OdorIndSmall,2,repnum(irep,5):repnum(irep,6)));        
    S = sum(squeeze(any(any(~isnan(A1),3),4)),3);
    pcutoff = (sp.p_cutoff./(sum(sp.OdorIndSmall)*S));
    pcutoff(S==0) = NaN;        
    A11 = squeeze(sum(sum(any(A1<repmat(pcutoff,[1 1 size(A1,3) size(A1,4) size(A1,5)]),4),3),5)); 
    indd = squeeze(sum(sum(any(~isnan(A1),4),3),5)); 
    sess = repmat([1:size(A11,2)],[size(A11,1) 1 size(A11,3)]);        
    ani = repmat([1:size(A11,1)],[size(A11,2) 1 size(A11,3)])';
    datA = A11./indd;        
    [p2,~,~] = anovan(datA(:),[sess(:) ani(:)],'varnames',{'Session';'Mouse'},'continuous',1,'display','off');
    
    if sum(ismember(figs_to_run,0))>0
        figure; hold on
        plot(sess(:),datA(:),'ko')
         if irep==1 %probes should decrease
            [r,p] = corr(sess(:),datA(:),'rows','complete','tail','left');
         elseif irep==2 %NTrepeats should increase
            [r,p] = corr(sess(:),datA(:),'rows','complete','tail','right');
         end
          N = sum(~isnan(sess(:)) & ~isnan(datA(:))); 
         suptitle(['Prop Corr Sig; N=' num2str(N) ' R = ' num2str(round(r,2,'significant')) ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' num2str(round(p2(1),3,'significant'))])         
         helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_PropCorrSig'])    
         
        A1 = squeeze(singledat.odorselective(:,:,:,sp.OdorIndSmall,2,repnum(irep,5):repnum(irep,6)));    
        A11 = squeeze(sum(any(A1<(sp.p_cutoff/sum(sp.OdorIndSmall)),4),3)); 
        indd = squeeze(sum(any(~isnan(A1),4),3)); 
        sess = repmat([1:size(A11,2)],[size(A11,1) 1 size(A11,3)]);        
        ani = repmat([1:size(A11,1)]',[1 size(A11,2) size(A11,3)]);
        repeatnum = permute(repmat([1:size(A11,3)]',[1 size(A11,1) size(A11,2)]),[2 3 1]);
        datA = A11./indd;        
        [p2,~,~] = anovan(datA(:),[sess(:) ani(:) repeatnum(:)],'varnames',{'Session';'Mouse';'RepeatNum'},'continuous',1,'display','off');
        
        figure; hold on
        plot(sess(:),datA(:),'ko')
         if irep==1 %probes should decrease
            [r,p] = corr(sess(:),datA(:),'rows','complete','tail','left');
         elseif irep==2 %NTrepeats should increase
            [r,p] = corr(sess(:),datA(:),'rows','complete','tail','right');
         end
          N = sum(~isnan(sess(:)) & ~isnan(datA(:))); 
         suptitle(['EachRepeatNum, N=' num2str(N) ' R = ' num2str(round(r,2,'significant')) ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' num2str(round(p2(1),3,'significant'))])         
         helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_EachRepeatNum'])    
         
        A1 = squeeze(singledat.odorselective(:,:,:,sp.OdorIndSmall,2,repnum(irep,5):repnum(irep,6)));    
        A11 = squeeze(sum(sum(any(A1<(sp.p_cutoff/sum(sp.OdorIndSmall)),4),3),5)); 
        indd = squeeze(sum(sum(any(~isnan(A1),4),3),5)); 
        sess = repmat([1:size(A11,2)],[size(A11,1) 1 size(A11,3)]);        
        ani = repmat([1:size(A11,1)]',[1 size(A11,2) size(A11,3)]);
        repeatnum = permute(repmat([1:size(A11,3)]',[1 size(A11,1) size(A11,2)]),[2 3 1]);
        datA = A11./indd;        
        [p2,~,~] = anovan(datA(:),[sess(:) ani(:) repeatnum(:)],'varnames',{'Session';'Mouse';'RepeatNum'},'continuous',1,'display','off');
        
        figure; hold on
        plot(sess(:),datA(:),'ko')
         if irep==1 %probes should decrease
            [r,p] = corr(sess(:),datA(:),'rows','complete','tail','left');
         elseif irep==2 %NTrepeats should increase
            [r,p] = corr(sess(:),datA(:),'rows','complete','tail','right');
         end
          N = sum(~isnan(sess(:)) & ~isnan(datA(:))); 
         suptitle(['AvgRepeat, N=' num2str(N) ' R = ' num2str(round(r,2,'significant')) ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' num2str(round(p2(1),3,'significant'))])         
         helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AvgRepeat'])    
    end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    if sum(ismember(figs_to_run,[0 16]))>0
        A1 = squeeze(singledat.odorselective(:,:,:,sp.OdorIndSmall,2,repnum(irep,1):repnum(irep,2)));
        A2 = squeeze(singledat.odorselective(:,:,:,sp.OdorIndSmall,2,repnum(irep,3):repnum(irep,4)));
        figa = figure; hold on
        figb = figure; hold on
        aa = NaN(size(A1,5),1); ylaa = aa;
        for ip = 1:size(A1,5)            
            inds = sum(sum(~isnan(A1(:,:,:,:,ip)),4)~=0,3);
            if sum(sum(sum(inds)))>0
                A11 = squeeze(sum(any(A1(:,:,:,:,ip)<(sp.p_cutoff/sum(sp.OdorIndSmall)),4),3)); 
                A22 = squeeze(sum(any(A2(:,:,:,:,ip)<(sp.p_cutoff/sum(sp.OdorIndSmall)),4),3)); 
                allp(:,:,1,irep,ip) = [A11./inds];
                allp(:,:,2,irep,ip) = [A22./inds];

                A11s = A11; A22s = A22;
                A11s(inds==0) = NaN; A22s(inds==0) = NaN;
                allp2(:,:,1,ip) = A11s;
                allp2(:,:,2,ip) = A22s;
                allp2(:,:,3,ip) = inds;

                figure(figa)
                aa(ip) = subplot(1,size(A1,5),ip); hold on
                for imouse = 1:4
                   plot(imouse,allp(imouse,:,1,irep,ip),'ok')                   
                   plot(imouse+6,allp(imouse,:,2,irep,ip),'or')
                end
                errorbar(2.5,mean(mean(allp(:,:,1,irep,ip),2,'omitnan'),1,'omitnan'),std(mean(allp(:,:,1,irep,ip),2,'omitnan'),1,'omitnan')./sqrt(size(allp,1)),'k')                
                errorbar(8.5,mean(mean(allp(:,:,2,irep,ip),2,'omitnan'),1,'omitnan'),std(mean(allp(:,:,2,irep,ip),2,'omitnan'),1,'omitnan')./sqrt(size(allp,1)),'r')
                set(gca,'xtick',[2.5 8.5],'xticklabel',{sp.rep_probe{irep},'Target'})      
                yl = get(gca,'ylim');
                ylaa(ip) = yl(2);
                xlim([0 11])
                ylabel('Proportion of neurons significant')                
                I1 = A11(:)./inds(:); I2 = A22(:)./inds(:);
                p = signrank(I1,I2,'tail','left');
                xlabel(['One-sided signrank: N =' num2str(sum(~isnan(I2))) ', P = ' num2str(round(p,2,'significant'))])

                figure(figb); subplot(1,size(A1,5),ip); hold on
                plot(1:2,[sum(A11,2) sum(A22,2)],'.-k','MarkerSize',20)  
                set(gca,'xtick',1:2,'xticklabel',{sp.rep_probe{irep},'Target'}) 
                xlim([.5 2.5])
                p2 = signrank(sum(A11,2),sum(A22,2),'tail','left');
                xlabel(['One-sided signrank: N =' num2str(size(A22,1)) ', P = ' num2str(round(p2,2,'significant'))])
            end
            title([sp.rep_probe{irep} ' ' num2str(ip) ' vs Target'])
        end

         figure(figa)
        for ip = 1:size(A1,5)            
            set(aa(ip),'ylim',[0 max(ylaa)])
%                 ylim([0 0.75])
        end
        
        A11 = squeeze(allp(:,:,1,irep,:));
        A22 = squeeze(allp(:,:,2,irep,:));
        
        p = signrank(A11(:),A22(:),'tail','left');
        suptitle(['One-sided signrank: N =' num2str(sum(~isnan(A11(:)))) ', P = ' num2str(round(p,2,'significant'))])
        set(gcf,'Position',[   1999         321         965         398])
        helper_saveandclosefig([thisdir 'PairedDownsample_' sp.ana_labs{3} '_' sp.rep_probe{irep}])

        % Figures for Liz's questions while editing the paper
        if sum(ismember(figs_to_run,0))>0
            A1 = squeeze(sum(allp2(:,:,1,:),2,'omitnan'));
            A2 = squeeze(sum(allp2(:,:,2,:),2,'omitnan'));
            A3 = squeeze(sum(allp2(:,:,3,:),2,'omitnan')); %total cells
    
         % Chi-squared test between up and down modulated and modulated by target at all
            % Observed data
            n1 = sum(A1(:),'omitnan'); N1 = sum(A3(:),'omitnan');
            n2 = sum(A2(:),'omitnan'); N2 = N1;
            % Pooled estimate of proportion
            p0 = (n1+n2) / (N1+N2);
            % Expected counts under H0 (null hypothesis)
            n10 = N1 * p0;
            n20 = N2 * p0;
            % Chi-square test, by hand
            observed = [n1 N1-n1 n2 N2-n2];
            expected = [n10 N1-n10 n20 N2-n20];
            [~,pp0,PPstats] = chi2gof([1 2 3 4],'freq',observed,'expected',expected,'ctrs',[1 2 3 4],'nparams',2);
    %         PPstats.Pvalue = pp0;
    
            figure(figb)
            p = signrank(A1(:),A2(:),'tail','left');
            suptitle(['One-sided signrank: N =' num2str(sum(~isnan(A11(:)))) ', P = ' num2str(round(p,2,'significant'))])
            set(gcf,'Position',[   1999         321         965         398])
            helper_saveandclosefig([thisdir 'PairedDownsample_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_forLiz'])
    
            figure; hold on;
            bar([observed(1:2)' observed(3:4)']','stacked')
            set(gca,'xtick',1:2,'xticklabel',{sp.rep_probe{irep},'Target'}) 
    %         xlim([.5 2.5])
            legend('Selective','Non-selective','Location','best')
            ylabel('Neurons')
            suptitle(['Chi-squared =' num2str(PPstats.chi2stat) ', P = ' num2str(round(pp0,2,'significant'))])
            helper_saveandclosefig([thisdir 'PairedDownsample_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_forLiz2'])
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if sum(ismember(figs_to_run,0))>0

        issig = any(any(singledat.odorsel(:,sp.OdorIndSmall,repnum(irep,5):repnum(irep,6),2)<(0.05/sum(sp.OdorIndSmall)),2),3);        
        isnotnan = any(any(~isnan(singledat.odorsel(:,sp.OdorIndSmall,repnum(irep,5):repnum(irep,6),2)),2),3);
        isnotnanodorsig = sd.all_Neurons(:,1) & any(any(~isnan(singledat.odorsel(:,sp.OdorIndSmall,repnum(irep,5):repnum(irep,6),2)),2),3);
        
        sess = NaN(size(sd.neuron_info,1),1);
        for isess = 1:max(sd.neuron_info(:,2))
            sess(sd.neuron_info(:,2)==isess,1) = sd.neuron_info(sd.neuron_info(:,2)==isess,3)-min(sd.neuron_info(sd.neuron_info(:,2)==isess,3))+1;
        end
        propsess = [];
        for isess = 1:max(sess)
            wrat = unique(sd.neuron_info(sess==isess,2));
            for irats = 1:length(wrat)
                ind = sess==isess & sd.neuron_info(:,2)==wrat(irats);
                propsess = [propsess; isess sum(issig(ind)) sum(isnotnan(ind)) sum(isnotnanodorsig(ind))];
            end
        end
        propsess2 = [propsess(propsess(:,3)>0,1) propsess(propsess(:,3)>0,2)./propsess(propsess(:,3)>0,3)];
         if irep==1 %probes should decrease
            [r,p] = corr(propsess2(:,1),propsess2(:,2),'rows','complete','tail','left');
         elseif irep==2 %NTrepeats should increase
            [r,p] = corr(propsess2(:,1),propsess2(:,2),'rows','complete','tail','right');
         end
         
        figure; hold on
        plot(propsess2(:,1),propsess2(:,2),'ko')
        N = sum(~isnan(propsess2(:,1)) & ~isnan(propsess2(:,2))); 
         suptitle(['N=' num2str(N) ' R = ' num2str(round(r,2,'significant')) ', P = ' num2str(round(p,3,'significant'))])
        ylabel(['Proportion of neurons significantly ' sp.rep_probe{irep} ' modulated'])
        xlabel(['Sessions'])
         helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance_FromStartSession'])    
         
          propsess2 = [propsess(propsess(:,4)>0,1) propsess(propsess(:,4)>0,2)./propsess(propsess(:,4)>0,4)];
         if irep==1 %probes should decrease
            [r,p] = corr(propsess2(:,1),propsess2(:,2),'rows','complete','tail','left');
         elseif irep==2 %NTrepeats should increase
            [r,p] = corr(propsess2(:,1),propsess2(:,2),'rows','complete','tail','right');
         end
         
        figure; hold on
        plot(propsess2(:,1),propsess2(:,2),'ko')
        N = sum(~isnan(propsess2(:,1)) & ~isnan(propsess2(:,2))); 
         suptitle(['N=' num2str(N) 'R = ' num2str(round(r,2,'significant')) ', P = ' num2str(round(p,3,'significant'))])
        ylabel(['Proportion of sig neurons significantly ' sp.rep_probe{irep} ' modulated'])
        xlabel(['Sessions'])
         helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance_FromStartSession_SigOdorOnly'])    
        
        %does the amount significant change across sessions?        
         inds = sum(sum(any(~isnan(A1(:,:,:,:,:))~=0,4),3),5);
         A11 = squeeze(sum(sum(any(A1(:,:,:,:,:)<(sp.p_cutoff/sum(sp.OdorIndSmall)),4),3),5)); 
         dat = A11./inds;
         dat2 = repmat(1:size(A11,2),[size(A11,1) 1]);
         figure; hold on;
         plot(dat2(:),dat(:),'ok')
         if irep==1 %probes should decrease
            [r,p] = corr(dat(:),dat2(:),'rows','complete','tail','left');
         elseif irep==2 %NTrepeats should increase
            [r,p] = corr(dat(:),dat2(:),'rows','complete','tail','right');
         end
         ylabel(['Proportion of neurons significantly ' sp.rep_probe{irep} ' modulated'])
         xlabel(['Sessions'])
         N = sum(~isnan(dat(:)) & ~isnan(dat2(:))); 
         suptitle(['N=' num2str(N) 'R = ' num2str(round(r,2,'significant')) ', P = ' num2str(round(p,3,'significant'))])
        
         helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_EachChance'])                
    end

    %figure 5d
    if (sum(ismember(figs_to_run,[0 5]))>0 && irep==2) || (sum(ismember(figs_to_run,0))>0)
         inds = sum(any(any(~isnan(A1(:,:,:,:,:))~=0,4),5),3);
         A11 = sum(any(any(A1(:,:,:,:,:)<(sp.p_cutoff/sum(sp.OdorIndSmall)),4),5),3);         
         dat = A11./inds;
         dat2 = repmat(1:size(A11,2),[size(A11,1) 1]);
       
        ani = repmat([1:size(A11,1)]',[1 size(A11,2) size(A11,3)]);
        [p2,tbl2,~] = anovan(dat(:),[dat2(:) ani(:)],'varnames',{'Session';'Mouse'},'continuous',1,'display','off');
    
        figure; hold on;
        for imouse = 1:4
            if sum(~isnan(dat(imouse,:)))>0
                p = plot(dat2(imouse,:),dat(imouse,:),'o-');
    %             [xF,yF] = helper_plotbestfitline(dat2(imouse,:),dat(imouse,:));
    %             plot(xF,yF,'Color',p.Color,'LineWidth',1)
            end
        end
        [xF,yF] = helper_plotbestfitline(dat2(:),dat(:));    
        plot(xF,yF,'Color','k','LineWidth',2)
    %     errorbar(mean(dat,1,'omitnan'),std(dat,[],1,'omitnan')./sqrt(sum(~isnan(dat),1)),'k')
        [r,p] = corr(dat(:),dat2(:),'rows','complete');
         ylabel(['Proportion of neurons significantly ' sp.rep_probe{irep} ' modulated'])
         xlabel(['Sessions'])
         N = sum(~isnan(dat(:)) & ~isnan(dat2(:))); 
         axis tight
         xl = get(gca,'xlim');
         set(gca,'xlim',[xl(1)-1 xl(2)+1])
         suptitle({['N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{irep}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})   
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance_mouse'])  
        writecell(tbl2,[thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance_mouse_ANOVAN.txt'],'Delimiter','tab')    
        writecell(tbl2,[thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance_mouse_ANOVAN.xlsx'])    
    end

    if sum(ismember(figs_to_run,0))>0

          ani = repmat([1:size(A11,1)]',[1 size(A11,2) size(A11,3)]);
        [p2,~,~] = anovan(dat(:),[dat2(:) ani(:)],'varnames',{'Session';'Mouse'},'continuous',1,'display','off');
        
    
    
        dat3 = NaN(size(dat2));
        [x,y] = find(~isnan(dat));
        for imouse = 1:max(x)
            dat3(imouse,y(x==imouse)) = 1:length(unique(y(x==imouse)));
        end
    
        [p2,~,~] = anovan(dat(:),[dat3(:) ani(:)],'varnames',{'Session';'Mouse'},'continuous',1,'display','off');
    
        figure; hold on;
        for imouse = 1:4
            if sum(~isnan(dat(imouse,:)))>0
                p = plot(dat3(imouse,:),dat(imouse,:),'o');
                [xF,yF] = helper_plotbestfitline(dat3(imouse,:),dat(imouse,:));
                plot(xF,yF,'Color',p.Color,'LineWidth',1)
            end
        end
        [xF,yF] = helper_plotbestfitline(dat3(:),dat(:));
        plot(xF,yF,'Color','k','LineWidth',2)
        [r,p] = corr(dat(:),dat3(:),'rows','complete');
         ylabel(['Proportion of neurons significantly ' sp.rep_probe{irep} ' modulated'])
         xlabel(['Sessions'])
         N = sum(~isnan(dat(:)) & ~isnan(dat3(:))); 
         suptitle({['N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{irep}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})   
        axis tight
              yl = get(gca,'ylim');
           xlim([0 19])       
    %      ylim([yl(1)-.1 yl(2)+.1])
        ylim([0 .7])
         set(gca,'xtick',5:5:15,'ytick',0:.2:1)
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance_mouse_SessionsInOrder'])         
    
    
        dat3 = NaN(size(dat2));
        [x,y] = find(~isnan(dat));
        for imouse = 1:max(x)
            dat3(imouse,1:sum(x==imouse)) = dat(imouse,~isnan(dat(imouse,:)));
            y(x==imouse) = 1:length(y(x==imouse));
        end
       
        figure; hold on;
        for imouse = 1:4
            if sum(~isnan(dat(imouse,:)))>0
                plot(dat3(imouse,:),'o-');
            end
        end
        ind = sum(sum(~isnan(dat3),1)>0);
        [xF,yF] = helper_plotbestfitline(1:ind,mean(dat3(:,1:ind),1,'omitnan'));
        plot(xF,yF,'Color','k','LineWidth',2)
        [r,p] = corr(y,dat3(~isnan(dat3)),'rows','complete');
         ylabel(['Proportion of neurons significantly ' sp.rep_probe{irep} ' modulated'])
         xlabel(['Sessions'])
         N = sum(~isnan(dat(:)) & ~isnan(dat3(:))); 
         suptitle({['N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{irep}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})   
        axis tight
              yl = get(gca,'ylim');
           xlim([0 19])       
    %      ylim([yl(1)-.1 yl(2)+.1])
        ylim([0 .7])
         set(gca,'xtick',5:5:15,'ytick',0:.2:1)
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance_mouse_SessionsInOrder2'])         
    
        
            figure; hold on;
         plot(dat2(:),dat(:),'ok') 
        [r,p] = corr(dat(:),dat2(:),'rows','complete');
    %      if irep==1 %probes should decrease
    %         [r,p] = corr(dat(:),dat2(:),'rows','complete','tail','left');
    %      elseif irep==2 %NTrepeats should increase
    %         [r,p] = corr(dat(:),dat2(:),'rows','complete','tail','right');
    %      end
         ylabel(['Proportion of neurons significantly ' sp.rep_probe{irep} ' modulated'])
         xlabel(['Sessions'])
         N = sum(~isnan(dat(:)) & ~isnan(dat2(:))); 
    %      suptitle(['N=' num2str(N) 'R = ' num2str(round(r,2,'significant')) ', P = ' num2str(round(p,3,'significant'))])
         suptitle({['N=' num2str(N) ', Sesions with ' ...
            sp.rep_probe{irep}];['R = ' num2str(round(r,2,'significant')) ...
            ' P = ' num2str(round(p,3,'significant')) ' ANOVAp = ' ...
            num2str(round(p2(1),3,'significant'))]})   
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance'])         
        
        inds = sum(sd.all_MeanSes(:,:,:,1)==1 & any(any(sum(~isnan(A1(:,:,:,:,:)),4)~=0,4),5),3);
        A11 = sum(sd.all_MeanSes(:,:,:,1)==1 & any(any(A1(:,:,:,:,:)<(sp.p_cutoff/sum(sp.OdorIndSmall)),4),5),3); 
        dat = A11./inds;
        figure; hold on
        plot(dat2(:),dat(:),'ok')
         if irep==1 %probes should decrease
            [r,p] = corr(dat(:),dat2(:),'rows','complete','tail','left');
         elseif irep==2 %NTrepeats should increase
            [r,p] = corr(dat(:),dat2(:),'rows','complete','tail','right');
         end
        ylabel(['Proportion of Sig neurons significantly ' sp.rep_probe{irep} ' modulated'])
        xlabel(['Sessions'])
         N = sum(~isnan(dat(:)) & ~isnan(dat2(:))); 
         suptitle(['N=' num2str(N) 'R = ' num2str(round(r,2,'significant')) ', P = ' num2str(round(p,3,'significant'))])
        helper_saveandclosefig([thisdir 'ProportionOverSessions_' sp.ana_labs{3} '_' sp.rep_probe{irep} '_AnyChance_SigOnly'])      
    end

end