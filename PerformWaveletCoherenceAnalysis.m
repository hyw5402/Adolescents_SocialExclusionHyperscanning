%% calculate the wavelet coherence Owen 20231023 
clear all; close all;
%% 
path_data='/home/DataDisk/Owen/fifth_hospital/Cyberball/AftPreprocess_CM_pca'; % the path containing NIRS data after preprocessing
path_cond='/home/DataDisk/Owen/fifth_hospital/Cyberball/Info';                 % the path contraining behavioral data
path_output='/home/DataDisk/Owen/fifth_hospital/Cyberball/WaveletCoherence_aftC_pca';  % the path for output
addpath('/home/DataDisk/Owen/fifth_hospital/Cyberball/Info');
if ~exist(path_output,'dir')
    mkdir(path_output);
end
save_originalRsq=0;

subjFiles=spm_select('List',path_data,'.mat');
subjFiles=cellstr(subjFiles); 




condition_data=load(fullfile(path_cond,'output_table_4wtc_2405.mat'));
condition_data=condition_data.output_table;

 parpool(30);  % open parallel workers

channel_s1=[1:16]; % channels selected; channel 17 & 18 are reference channels
channel_s2=[19:34];
num_fs=80;

 parfor i=1:length(subjFiles)
 % for i=1%:length(subjFiles)   
 % read data file
    nirs_data=load(fullfile(path_data,subjFiles{i}));
    % nirs_data=nirs_data.nirsdata;
    nirs_data=nirs_data.data_filter;
    % read condition file
    temp_name=strrep(subjFiles{i},'.mat','');
    temp_cond_data=condition_data(temp_name,:);
    temp_cond_data=table2array(temp_cond_data);
    
%     mean_indv_wlc_Rsq=zeros(length(combined_1),length(combined_2),5);
    mean_indv_wlc_Rsq=cell(length(channel_s1),length(channel_s2),5);  % initiate output cells
    
    mean_indv_wlc_Rsq_3stage=cell(length(channel_s1),length(channel_s2),3);
    
    % indv_wlc_Rsq=cell(length(channel_s1),length(channel_s2),5);
    % indv_wlc_coi=cell(length(channel_s1),length(channel_s2),5);
    indv_wlc_period=cell(1);
    error_status=zeros(length(channel_s1),length(channel_s2));
    error_messages=cell(length(channel_s1),length(channel_s2));
%     indv_wlc.coi=cell(length(combined_1),length(combined_2),5);
%     indv_wlc.sig95=cell(length(combined_1),length(combined_2),5);
    task_onset=ceil(temp_cond_data(1)*10+15*10);
    task_end=floor(temp_cond_data(10)*10);
    if task_end>length(nirs_data)
        task_end=length(nirs_data);
    end
        
        n=0;
        for j=1:length(channel_s1)
            for jj=1:length(channel_s2)
            
%        
            x_data_1=nirs_data(task_onset:task_end,channel_s1(j));
            y_data_1=nirs_data(task_onset:task_end,channel_s2(jj));
            
             
            % wlc

                   try
          
             % [Rsq_1,period_1,~,coi_1,~]=wtc(x_data_1,y_data_1,'S0',10,'ms',1000,'AR1','auto'); 
             % incoi_1=period_1(:)*(1./coi_1)>1;
                  [Rsq_1,~,period_1,coi_1,~]=wcoherence(x_data_1,y_data_1,10,'FrequencyLimits',[0.01,1]); % calculate wavelet coherence in frequency range: 0.01-1
            
                  % reset abnormal data
                   Rsq_1(Rsq_1<0)=0;
                   Rsq_1(Rsq_1>=1)=0.99;

              % fisher z
                  Rsq_1= 0.5.*log((1+Rsq_1)./(1-Rsq_1));
                  
% mean_indv_wlc_Rsq(j,jj,k)=mean(Rsq,'all','omitnan');

                
                
                if j==1 && jj==1
            indv_wlc_period{1}=period_1;
          
                end
              catch error_msg
    error_status(j,jj)=1;
    error_messages{j,jj}=error_msg;
              end
            % incoi=period(:)*(1./coi)>1;
       % 

       %% arrange results
    for k=1:5
      
        if k==1
            onset_1=1;
        else
        onset_1=temp_cond_data(2*k-1)*10-temp_cond_data(1)*10-15*10+1;
        end
      
       onset_2=temp_cond_data(2*k)*10-temp_cond_data(1)*10-15*10;
       
       if onset_2>size(Rsq_1,2)
           onset_2=size(Rsq_1,2);
       end
     temp_rsq=mean(Rsq_1(:,onset_1:onset_2),2,'omitnan');
  mean_indv_wlc_Rsq{j,jj,k}=temp_rsq;

    if k<4
        if k==1
         onset_1_3stage=temp_cond_data(2*k-1)*10-temp_cond_data(1)*10+1;
       
      
       onset_2_3stage=temp_cond_data(2*k)*10-temp_cond_data(1)*10-15*10;
        elseif k==2
          onset_1_3stage=temp_cond_data(2*k-1)*10-temp_cond_data(1)*10-15*10+1;
       
      
       onset_2_3stage=temp_cond_data(2*k)*10-temp_cond_data(1)*10-15*10;
        elseif k==3
            onset_1_3stage=temp_cond_data(2*3-1)*10-temp_cond_data(1)*10-15*10+1;
       
      
       onset_2_3stage=temp_cond_data(2*5)*10-temp_cond_data(1)*10-15*10;
      if onset_2_3stage>size(Rsq_1,2)
           onset_2_3stage=size(Rsq_1,2);
       end

        end
        temp_rsq_3stage=mean(Rsq_1(:,onset_1_3stage:onset_2_3stage),2,'omitnan');
  mean_indv_wlc_Rsq_3stage{j,jj,k}=temp_rsq_3stage;
    end
       
         
    end
         

            end
            
        end
     



       % fprintf('\n %s condition %d finished!',temp_name);
%     output_1{i}=mean_indv_wlc_Rsq;
%     output_2{i}=indv_wlc_period;
      output_data=cell(5,2);
      output_data(:,1)={'period';'mean_indv_wlc_Rsq';'mean_indv_wlc_Rsq_3stage';'error_status';'error_messages'};
%        if save_originalRsq==1
%       output_data{1,2}=indv_wlc_Rsq;
% %       output_data{2,2}=indv_wlc_coi;
%        end
      output_data{1,2}=indv_wlc_period;
      output_data{2,2}=mean_indv_wlc_Rsq;
      output_data{3,2}=mean_indv_wlc_Rsq_3stage;

      output_data{4,2}=error_status;
      output_data{5,2}=error_messages;
%       temp_name=strrep(subjFiles{i},'_converted_data.mat','');
      parsave(fullfile(path_output,[temp_name '.mat']),output_data); % save data
    fprintf('\n %s finished!',subjFiles{i});
end
        

    