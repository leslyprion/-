%%
% 
% <<LaserModulatedUnit_summary_Raster_PSTH.tif>>
% 
% 
% 

% clear;close all;
[all_mat, FILE_PATH] = uigetfile('*.mat', 'Choose MAT File', 'multiSelect', 'on');
if ~iscell(all_mat)
    all_mat = {all_mat}; 
end
cd(FILE_PATH);
nUnit = size(all_mat, 2);
Plot_baseline = -1;
Plot_xlim = 1.5;
Trial_sort_start = 150;
Trial_sort_end = 220;
fig_append = ['Trial',num2str(Trial_sort_start),'to',num2str(Trial_sort_end)];

Laser_sig = 0;
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 2 42 22],'Visible','on');
for iUnit = 1:nUnit
    load(all_mat{1,iUnit});
    chan_unit = [vertcat(Unit_FRchangedbylaser_AUC.Chan),vertcat(Unit_FRchangedbylaser_AUC.Unit)];
    chan_unit_index = find(chan_unit(:,1)==UnitChanInfo(1,2) & chan_unit(:,2) == UnitChanInfo(1,3));
    if ~isempty(chan_unit_index)
        Laser_sig = Laser_sig+1;
        Stimulus = NeuralActiByTrial(Trial_sort_start:Trial_sort_end);
        Trial = TrialInformation(Trial_sort_start:Trial_sort_end);
        CR_LaserOn_trial = find(vertcat(Trial.TrialOutput)==3 & vertcat(Trial.Laser)==1);
        CR_LaserOff_trial = find(vertcat(Trial.TrialOutput)==3 & vertcat(Trial.Laser)==0);
        Hit_LaserOn_trial = find(vertcat(Trial.TrialOutput)==1 & vertcat(Trial.Laser)==1);
        Hit_LaserOff_trial = find(vertcat(Trial.TrialOutput)==1 & vertcat(Trial.Laser)==0);
        Miss_LaserOn_trial = find(vertcat(Trial.TrialOutput)==5 & vertcat(Trial.Laser)==1);
        Miss_LaserOff_trial = find(vertcat(Trial.TrialOutput)==5 & vertcat(Trial.Laser)==0);
        FA_LaserOn_trial = find(vertcat(Trial.TrialOutput)==4 & vertcat(Trial.Laser)==1);
        FA_LaserOff_trial = find(vertcat(Trial.TrialOutput)==4 & vertcat(Trial.Laser)==0);
        Nogo_LaserOn_trial = [CR_LaserOn_trial;FA_LaserOn_trial];
        Go_LaserOn_trial = [Hit_LaserOn_trial;Miss_LaserOn_trial];
        Nogo_LaserOff_trial = [CR_LaserOff_trial;FA_LaserOff_trial];
        Go_LaserOff_trial = [Hit_LaserOff_trial;Miss_LaserOff_trial];
        
        Go_LaserOff_spike = Stimulus(Go_LaserOff_trial);
        Go_LaserOn_spike = Stimulus(Go_LaserOn_trial);
        Nogo_LaserOff_spike = Stimulus(Nogo_LaserOff_trial);
        Nogo_LaserOn_spike = Stimulus(Nogo_LaserOn_trial);
        
        LaserOn_GNG_Spike = [Stimulus(Go_LaserOn_trial),Stimulus(Nogo_LaserOn_trial)];
        LaserOff_GNG_Spike = [Stimulus(Go_LaserOff_trial),Stimulus(Nogo_LaserOff_trial)];
        
        
        SAMPLING_RATE = Info_nsx.Fs_Hz;
        BIN_SIZE = BasicER.PSTHbin;
        TimeStamp = -TrialInformation(1).Baseline/1000+(BIN_SIZE/SAMPLING_RATE/2):BIN_SIZE/SAMPLING_RATE:(Basic.Delay_duration+Basic.Stimulus_duration)/1000-(BIN_SIZE/SAMPLING_RATE/2);

        %% 
%         fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 2 42 22],'Visible','on');
        %%%%%%%%%%%%%%%%%%%%%
        m = ceil(Laser_sig/4);
        n = mod(Laser_sig,4);if n == 0; n = 4;end
        axOff_Raster = subplot(6,8,(m-1)*16+(n-1)*2+1);hold(axOff_Raster,'on');
        index_spike = Go_LaserOff_spike;
        for iTrial = 1:length(index_spike)
            if ~isempty(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)
            index_raster = errorbar_2(index_spike(iTrial).BaseCueDelaySpikeTimeRelative/SAMPLING_RATE,iTrial*ones(size(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)),...
                0.7*ones(size(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)));
            set(index_raster(2),'linestyle','none');
            set(index_raster(1),'color','b');
            end
        end
        y_max = length(index_spike);
        index_spike = Nogo_LaserOff_spike;
        for iTrial = 1:length(index_spike)
            if ~isempty(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)
            index_raster = errorbar_2(index_spike(iTrial).BaseCueDelaySpikeTimeRelative/SAMPLING_RATE,(iTrial+5+y_max)*ones(size(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)),...
                0.7*ones(size(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)));
            set(index_raster(2),'linestyle','none');
            set(index_raster(1),'color','r');
            end
        end
        set(axOff_Raster,'XLim',[Plot_baseline Plot_xlim]);
        title(axOff_Raster,{['Chan',num2str(UnitChanInfo(1,2)),'Unit',num2str(UnitChanInfo(1,3))];'Laser Off'});
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        axOff_PSTH = subplot(6,8,(m-1)*16+(n-1)*2+1+8);hold(axOff_PSTH,'on');
        index_PSTH = vertcat(Go_LaserOff_spike.BaseCueDelaySpikePsth);
        index_PSTH_mean = mean(index_PSTH,1);
        index_errorbar = std(index_PSTH,1)/sqrt(size(index_PSTH,1));
        plot(axOff_PSTH,TimeStamp,index_PSTH_mean,'-','Color','b','DisplayName','Go');
        F1 = fill(axOff_PSTH,[TimeStamp fliplr(TimeStamp)],[index_PSTH_mean-index_errorbar fliplr(index_PSTH_mean+index_errorbar)],'b','FaceAlpha',0.3,'EdgeAlpha',0);
        set(get(get(F1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        index_PSTH = vertcat(Nogo_LaserOff_spike.BaseCueDelaySpikePsth);
        index_PSTH_mean = mean(index_PSTH,1);
        index_errorbar = std(index_PSTH,1)/sqrt(size(index_PSTH,1));
        plot(axOff_PSTH,TimeStamp,index_PSTH_mean,'-','Color','r','DisplayName','Nogo');
        F1 = fill(axOff_PSTH,[TimeStamp fliplr(TimeStamp)],[index_PSTH_mean-index_errorbar fliplr(index_PSTH_mean+index_errorbar)],'r','FaceAlpha',0.3,'EdgeAlpha',0);
        set(get(get(F1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        if Laser_sig == 1
         L = legend(axOff_PSTH);set(L,'box','off','Location','NorthWest');
        end
        set(axOff_PSTH,'XLim',[Plot_baseline Plot_xlim]);
%         xlabel(axOff_PSTH,'Time(s)');ylabel(axOff_PSTH,'FR (Hz)');
        yaxis1 = get(axOff_PSTH,'YLim');
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        axOn_Raster = subplot(6,8,(m-1)*16+(n-1)*2+1+1);hold(axOn_Raster,'on');
        index_spike = Go_LaserOn_spike;
        for iTrial = 1:length(index_spike)
            if ~isempty(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)
            index_raster = errorbar_2(index_spike(iTrial).BaseCueDelaySpikeTimeRelative/SAMPLING_RATE,iTrial*ones(size(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)),...
                0.7*ones(size(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)));
            set(index_raster(2),'linestyle','none');
            set(index_raster(1),'color',[0.07,0.62,1.00]);
            end
        end
        y_max = length(index_spike);
        index_spike = Nogo_LaserOn_spike;
        for iTrial = 1:length(index_spike)
            if ~isempty(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)
            index_raster = errorbar_2(index_spike(iTrial).BaseCueDelaySpikeTimeRelative/SAMPLING_RATE,(iTrial+5+y_max)*ones(size(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)),...
                0.7*ones(size(index_spike(iTrial).BaseCueDelaySpikeTimeRelative)));
            set(index_raster(2),'linestyle','none');
            set(index_raster(1),'color',[1.00,0.41,0.16]);
            end
        end
        set(axOn_Raster,'XLim',[Plot_baseline Plot_xlim]);
        title(axOn_Raster,{['Chan',num2str(UnitChanInfo(1,2)),'Unit',num2str(UnitChanInfo(1,3))];'Laser On'});
        
        %%%%%%%%%%%%%%%%%%%%%
        axOn_PSTH = subplot(6,8,(m-1)*16+(n-1)*2+1+1+8);hold(axOn_PSTH,'on');
        index_PSTH = vertcat(Go_LaserOn_spike.BaseCueDelaySpikePsth);
        index_PSTH_mean = mean(index_PSTH,1);
        index_errorbar = std(index_PSTH,1)/sqrt(size(index_PSTH,1));
        plot(axOn_PSTH,TimeStamp,index_PSTH_mean,'-','Color',[0.07,0.62,1.00],'DisplayName','Go');
        F1 = fill(axOn_PSTH,[TimeStamp fliplr(TimeStamp)],[index_PSTH_mean-index_errorbar fliplr(index_PSTH_mean+index_errorbar)],[0.07,0.62,1.00],'FaceAlpha',0.3,'EdgeAlpha',0);
        set(get(get(F1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        index_PSTH = vertcat(Nogo_LaserOn_spike.BaseCueDelaySpikePsth);
        index_PSTH_mean = mean(index_PSTH,1);
        index_errorbar = std(index_PSTH,1)/sqrt(size(index_PSTH,1));
        plot(axOn_PSTH,TimeStamp,index_PSTH_mean,'-','Color',[1.00,0.41,0.16],'DisplayName','Nogo');
        F1 = fill(axOn_PSTH,[TimeStamp fliplr(TimeStamp)],[index_PSTH_mean-index_errorbar fliplr(index_PSTH_mean+index_errorbar)],[1.00,0.41,0.16],'FaceAlpha',0.3,'EdgeAlpha',0);
        set(get(get(F1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%         L = legend(axOn_PSTH);set(L,'box','off','Location','NorthWest');
        set(axOn_PSTH,'XLim',[Plot_baseline Plot_xlim]);
%         xlabel(axOn_PSTH,'Time(s)');ylabel(axOn_PSTH,'FR (Hz)');
        yaxis2 = get(axOn_PSTH,'YLim');
        
        yaxis = [min(yaxis1(1),yaxis2(1)) max(yaxis1(2),yaxis2(2))];
        set(axOn_PSTH,'YLim',yaxis);
        set(axOff_PSTH,'YLim',yaxis);
        
        sgtitle(fig_append);
        
        
        
        
    
    end
end