function singledat = calculate_plot_singleneuron_selectivity(dirs, ...
    neuron_info,neuronlistnew,mousenames,toplotcells,lookat, ...
    labeladd,odorwindowlabel,sm,num_avg_over)

%%%%%
%%%%%
%%%%% This script is to calculate for each neuron, recorded across many
%%%%% sessions with different conditions and trial types, what, if
%%%%% anything, it is significantly modulated to. This includes odor onset,
%%%%% odor offset, licks, and also individual or groups of odors (target,
%%%%% non-target, probe, non-target repeat). Then group analyses are
%%%%% performed to see the proportions of these categories (e.g. target
%%%%% selective neurons), their temporal dynamics in the trial, and how the
%%%%% proportion of this type of neuron changes across sessions. 
%%%%%
%%%%% It also includes determining whether neurons are significantly
%%%%% Target/NT modulated when the number of trials in each group is 
%%%%% down-sampled to 2/3ed of probe or NTrepeat trials depending on
%%%%% which they are being compared to. This 2/3rds is chosen randomly 100
%%%%% times and the mode of the p-values generated by those 100 comparisons
%%%%% is chosen as the p-value. This part adds a lot of time. It
%%%%% should take only ~20min to generate all the example neurons needed in
%%%%% the paper, but to run the orginial significance numbers for all the
%%%%% cells it takes much longer.
%%%%%
%%%%%

%{For each trial the preodor is 5 seconds, odor on is 2 seconds, and 
% post odor is 10 seconds. (i.e. raster_data is 290 x 17000 where each 
% row is 1 of 290 trials with 17000 ms trial duration) get rewarded for 
% licks after solenoid opens (0ms in lick data, 5000 in raster 
% data... so ~5500 in raster data is when rewarded licks start happening)}

try

%%% set up parameters 
tt_all = 0;
Ls = [111 222 333; 11 22 33]; % probe and nontarget markers
numbins = size(odorwindowlabel,1);
odorwindow = odorwindowlabel+5000;
odorwindowL = odorwindowlabel;
% first 4 seconds and last 4 seconds is baseline
baseline = [0 4000; 13000 17000]; 
numshuff = 1000;
% to isolate trials that had a large difference between the lick and the 
% odor-on        
sepcut = 800; 
if isempty(lookat)
    lookat = 1:size(neuron_info,1);
    %debug lookat = 816:size(neuron_info,1);
end

%%% set up figure labels
figlablab = {'Target';'Probe1';'Probe2';'Probe3';'Repeat1';'Repeat2'};
labeladd = [labeladd '_numshuff' num2str(numshuff)];

%%% set up variables to fill
sessionneurons = zeros(4,38); 
 %rats, days, neurons, bin, modulation/p-value/rawavgLicktrials/rawavg,
 % odor/lick/target/correct/probe/repeat / compare target with fewer trials
odorselective = NaN(4,38,30,numbins,4,19);
 %odor/lick/target/correct/probex3/repeatx2
odorsel = NaN(size(neuron_info,1),numbins,19,4);
%mean/median clasify, probe/repeat, #

 %target/correct/off response target/ 
 % off response correct / overall odor-on / overall odor-off
odorsel_all = NaN(size(neuron_info,1),6,3);
odorselective_all = NaN(4,38,30,6,3);
singledat.toplotOdorLick = NaN(size(neuron_info,1),12000);
neuron_info(:,6) = NaN(size(neuron_info,1),1);


%%% Get data out and plot individual example neurons
for ineuron2 = 1:length(lookat) 
    
    tic

    ineuron = lookat(ineuron2);
    
    % get session data and what neuron of the session and add to neuron_info
    sessionneurons(neuron_info(ineuron,2),neuron_info(ineuron,3)) = ...
        sessionneurons(neuron_info(ineuron,2),neuron_info(ineuron,3)) +1;
    neuronnum =  ...
        sessionneurons(neuron_info(ineuron,2),neuron_info(ineuron,3));
    neuron_info(ineuron,6) = neuronnum;
    
    % load data
    load(neuronlistnew{ineuron},'raster_labels','raster_data')
  
    % get trial types and correct/incorrect
    [trialtypes, cor, isProbe, isRepeat] = get_trial_info(raster_labels);
    
    %%% get neural data

    % odor-on windows
     %500ms before odor-on
    predat = sum(raster_data(:,4501:5000),2)./.5;
     %100-600ms after odor-on
    postdat = sum(raster_data(:,5101:5600),2)./.5; 
    
    % odor-off windows
    predat_off = sum(raster_data(:,6501:7000),2)./.5;
    postdat_off = sum(raster_data(:,7101:7600),2)./.5; 
    %100-600ms after odor-off
        
    % make target/nt only
    reallabelsTarget = NaN(size(trialtypes))';
    reallabelsTarget(trialtypes==1) = 1;
    reallabelsTarget(trialtypes==2) = 2;
    
    % overall pvalues of target, correct/error for odor-on windows   
    [Target_p, Target_r] = ...
        get_shuffle_p_single_neuron_diff(postdat-predat,...
        reallabelsTarget,numshuff,1:2);
    [Correct_p, Correct_r] = ...
        get_shuffle_p_single_neuron_diff(postdat-predat,...
        cor',numshuff,[1 0]);  

    % p-values overall for odor-off periods
    [Target_p_off, Target_r_off] = ...
        get_shuffle_p_single_neuron_diff(postdat_off-predat_off,...
        reallabelsTarget,numshuff,1:2);
    [Correct_p_off, Correct_r_off] = ...
        get_shuffle_p_single_neuron_diff(postdat_off-predat_off,...
        cor',numshuff,[1 0]);

                        
    %get data in different format
    [x,y] = find(raster_data);       
    ysh = NaN(length(y),numshuff);
    for iE = 1:size(raster_data,1)    
        k = repmat((randi(size(raster_data,2),numshuff,1))',...
            [sum(x==iE) 1]); 
        j = repmat(y(x==iE),[1 numshuff]);
        j2 = mod(j+k,size(raster_data,2));
        ysh(x==iE,:) = j2;
    end

    %Licks
    thelicktime = get_licks(raster_labels);

    % get index to make the lick on each trial at 0
    yL = NaN(size(y)); yL2 = NaN(size(y,1),3);
    yLsh = NaN(length(y),numshuff);
    for ix = 1:max(x)
        yL(x==ix) = y(x==ix)-(thelicktime(ix)+5000);
        yLsh(x==ix,:) = ysh(x==ix,:)-(thelicktime(ix)+5000);
        yL2(x==ix,:) = [y(x==ix) y(x==ix)-thelicktime(ix) ...
            thelicktime(ix)*ones(sum(x==ix),1)];
    end    
    
    %overall p-value using the normalized modulation of odor-on and
    %odor-off periods
    [p_overall, mod1_overall] = get_shuffle_p_single_neuron_modu(y,yL,...
        [4501 5000],[5101 5600],ysh,[]);
    [p_overall_off, mod1_off] = get_shuffle_p_single_neuron_modu(y,yL,...
        [6501 7000],[7101 7600],ysh,[]);

    
    %get p-values and modulation for each bin post odor onset and odor
    %offset
    for ibin = 1:numbins
        
        [p1, mod1, odorbin] = get_shuffle_p_single_neuron_modu(y,yL,...
            baseline,odorwindow(ibin,:),ysh,[]);
        [p2, modlick1, lickbin] = get_shuffle_p_single_neuron_modu(y,yL,...
            baseline,odorwindowL(ibin,:),ysh,yLsh);     

        allmod = mean(raster_data(:,...
            odorwindow(ibin,1)+1:odorwindow(ibin,2)),2,'omitnan');                       
        [Target_pBin, Target_rBin] = get_shuffle_p_single_neuron_diff...
            (allmod,reallabelsTarget,numshuff,1:2);
        [Corr_pBin, Corr_rBin] = get_shuffle_p_single_neuron_diff...
            (allmod,cor',numshuff,[1 0]);

        
        Probe_rBin = NaN(1,3); Probe_pBin2 = Probe_rBin;
        Probe_pBin = NaN(1,3);        
        Target_pBin_compareProbe = NaN(1,3);
        if isProbe
            for ip = 1:3
                if sum(trialtypes==Ls(1,ip))>0

                   [Probe_pBin(ip),Probe_rBin(ip)] = ...
                       get_shuffle_p_single_neuron_diff...
                       (allmod,trialtypes,numshuff,[Ls(1,ip) 2]);
                   
                   % statistics down-sampling to the same number (shuffled
                   % num_avg_over times and mode is taken)
                   num_reduce_probe = floor(sum(trialtypes==Ls(1,ip))./1.5);
                   
                   Probe_pBin_temp = NaN(num_avg_over,1);                   
                   testind = find(trialtypes==Ls(1,ip));
                   for ipp = 1:num_avg_over
                       randtestind = trialtypes;
                       randtestind(testind(randperm(length(testind),...
                           length(testind)-num_reduce_probe))) = NaN;  
                       Probe_pBin_temp(ipp) = ...
                           get_shuffle_p_single_neuron_diff(allmod,...
                           randtestind,numshuff,[Ls(1,ip) 2]);
                   end
                   Probe_pBin2(ip) = mode(round(Probe_pBin_temp*1000)/1000);
                   
                    testind = find(reallabelsTarget==1);                  
                    Target_pBin_compareProbe_temp = NaN(num_avg_over,1);
                    for ipp = 1:num_avg_over
                        randtestind = trialtypes; 
                        randtestind(testind(randperm(length(testind),...
                            length(testind)-num_reduce_probe))) = NaN;

                        Target_pBin_compareProbe_temp(ipp) = ...
                            get_shuffle_p_single_neuron_diff(allmod,...
                            randtestind,numshuff,[1 2]);
                    end
                    Target_pBin_compareProbe(ip) = ...
                        mode(round(Target_pBin_compareProbe_temp*1000)/1000);
                    
                else
                    Probe_pBin(ip) = NaN; 
                    Probe_rBin(ip) = NaN; 
                    Probe_pBin2(ip) = NaN;
                end
            end
        end
        
        Repeat_pBin = NaN(1,2);
        Repeat_rBin = NaN(1,2);     Repeat_pBin2 = Repeat_pBin;   
        Target_pBin_compareRepeat = NaN(1,2);
        if isRepeat         
            for ip = 1:2
                if sum(trialtypes==Ls(2,ip))>0
                                        
                   [Repeat_pBin(ip), Repeat_rBin(ip)] = ...
                       get_shuffle_p_single_neuron_diff...
                       (allmod,trialtypes,numshuff,[Ls(2,ip) 2]);
                    
                   % statistics down-sampling to the same number (shuffled
                   % num_avg_over times and averaged)
                   num_reduce_repeat = floor(sum(trialtypes==Ls(2,ip))./1.5);                    
                   Repeat_pBin_temp = NaN(num_avg_over,1);                   
                   testind = find(trialtypes==Ls(2,ip));
                   for ipp = 1:num_avg_over
                       randtestind = trialtypes;
                       randtestind(testind(randperm...
                           (length(testind),length(testind)-...
                           num_reduce_repeat))) = NaN;     

                       Repeat_pBin_temp(ipp) = ...
                           get_shuffle_p_single_neuron_diff...
                           (allmod,randtestind,numshuff,[Ls(2,ip) 2]);
                   end
                   Repeat_pBin2(ip) = ...
                       mode(round(Repeat_pBin_temp*1000)/1000);
                                   
                   testind = find(reallabelsTarget==1);      
                   Target_pBin_compareRepeat_temp = NaN(num_avg_over,1);
                   for ipp = 1:num_avg_over
                       randtestind = trialtypes; 
                       randtestind(testind(randperm...
                           (length(testind),length(testind)-...
                           num_reduce_repeat))) = NaN;        

                       Target_pBin_compareRepeat_temp(ipp) = ...
                           get_shuffle_p_single_neuron_diff...
                           (allmod,randtestind,numshuff,[1 2]);
                   end
                   Target_pBin_compareRepeat(ip) = ...
                       mode(round(Target_pBin_compareRepeat_temp*1000)/1000);

                end
            end
        end         
        
        % to isolate trials that had a large difference between the 
        % lick and the odor-on use sepcut
        odorbin2 = (sum(y(~isnan(yL) & yL2(:,3) > sepcut) > ...
            odorwindow(ibin,1) & y(~isnan(yL) & yL2(:,3)>sepcut) <= ...
            odorwindow(ibin,2))./diff(odorwindow(ibin,:)./1000));
%             ./sum(~isnan(thelicktime) & thelicktime>sepcut); 
        %length(unique(x)); %taking the FR for all trials with at 
        % least one spike (could change)

        lickbin2 = (sum(yL>odorwindowL(ibin,1) & yL<=odorwindowL(ibin,2)...
            & yL2(:,3)>sepcut)./diff(odorwindowL(ibin,:)./1000)); 
        %./sum(~isnan(thelicktime) & thelicktime>sepcut); 
        %length(unique(x));
       
        % for plotting
        odorsel(ineuron,ibin,:,1) = [mod1 modlick1 Target_rBin Corr_rBin...
            Probe_rBin Repeat_rBin NaN(1,5) Probe_rBin Repeat_rBin];
        odorsel(ineuron,ibin,:,2) = [p1 p2 Target_pBin Corr_pBin ...
            Probe_pBin2 Repeat_pBin2 Target_pBin_compareProbe ...
            Target_pBin_compareRepeat Probe_pBin Repeat_pBin];
        odorsel(ineuron,ibin,:,3) = [odorbin2 lickbin2 NaN(1,17)];
        odorsel(ineuron,ibin,:,4) = [odorbin lickbin NaN(1,17)];
        
        % make this into a structure %%%%%%
        %rats, days, neurons, bin, p-value/modulation/FR/FRin large lick and odor diff trials, odor/lick/target/correct                
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,1) = [mod1 p1 odorbin2 odorbin]; 
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,2) = [modlick1 p2 lickbin2 lickbin];
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,3) = [Target_rBin Target_pBin NaN(1,2)];      
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,4) = [Corr_rBin Corr_pBin NaN(1,2)];
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,5:7) = [Probe_rBin; Probe_pBin2; NaN(2,3)];             
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,8:9) = [Repeat_rBin; Repeat_pBin2; NaN(2,2)];     
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,10:12) = [NaN(1,3); Target_pBin_compareProbe; NaN(2,3)];
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,13:14) = [NaN(1,2); Target_pBin_compareRepeat; NaN(2,2)];        
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,15:17) = [Probe_rBin; Probe_pBin; NaN(2,3)];             
        odorselective(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,ibin,:,18:19) = [Repeat_rBin; Repeat_pBin; NaN(2,2)];  
    end
    
    odorsel_all(ineuron,:,1) = [Target_r Correct_r Target_r_off ...
        Correct_r_off mod1_overall mod1_off];
    odorsel_all(ineuron,:,2) = [Target_p Correct_p Target_p_off ...
        Correct_p_off p_overall p_overall_off];
    
    odorselective_all(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,:,1) = [Target_r Correct_r Target_r_off Correct_r_off mod1_overall mod1_off];
    odorselective_all(neuron_info(ineuron,2),neuron_info(ineuron,3),neuronnum,:,2) = [Target_p Correct_p Target_p_off Correct_p_off p_overall p_overall_off];
   
    % Get time data to plot either here or with group
    [x1,~] = find(raster_data);  
    [~,ord] = sort(thelicktime);
    
    x = NaN(size(x1));
    for ix = 1:max(x1)
        x(x1==ord(ix)) = ix;            
    end

    %%%% what is this for
    toplotOdor = raster_data(:,4000:10000);
    toplotLick = zeros(size(raster_data));
    for ix = 1:length(x1)
        if yL(ix)>-6000 && yL(ix)<10000
            toplotLick(x1(ix),yL(ix)+6000) = 1;  
        end
    end
    toplotLick(isnan(thelicktime),:) = [];    
    singledat.toplotOdorLick(ineuron,:) = cat(2,smoothdata(mean(toplotOdor(:,1:end-1)),'movmean',sm),...
        smoothdata(mean(toplotLick(:,4000:10000-1)),'movmean',sm));    

    if toplotcells
        %could streamline this more
        plot_singleneuron_selectivity(Ls,ineuron,yL,raster_data,dirs,...
            odorsel,sm,toplotOdor,figlablab,mousenames,labeladd,...
            odorwindowlabel,toplotLick,neuron_info,neuronnum,...
            thelicktime,reallabelsTarget,cor,odorsel_all,...
            isProbe,isRepeat,trialtypes)
    end

%     if rem(ineuron,10)==0
    tt = toc;
    tt_all = tt_all+(tt/60);
    disp(['neuron ' num2str(ineuron2) ' of ' num2str(length(lookat)) ' (' num2str(round(100*ineuron2/length(lookat),1,"decimals")) '%) in ' num2str(round(tt/60,2,'significant')) ' min (' num2str(round(tt_all,2,'significant')) ' total)'    ])
%     end
        
end

singledat.odorsel_all = odorsel_all;
singledat.odorselective_all = odorselective_all;
singledat.neuron_info = neuron_info;
singledat.odorsel = odorsel;
singledat.odorselective = odorselective;


catch ME
   ME.getReport
   disp('*')
end