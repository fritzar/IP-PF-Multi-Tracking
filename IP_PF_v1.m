% ==================
clear all;
close all;
clc;
format long; 

%% ======The initialization of the basic parameters=====
Num_Cell_x=50;
Num_Cell_y=50;

Total_time=37;     % the number of total integrated frames
Total_time_data=37;
Time_point=37;
T_step=1;          % The size of the time cell:Time_step

q1=0.001;          % q1,q2 the level of process noise in target motion
q1_1=0.01;

Re_x = 1; %Resolution of x-axis
Re_y = 1;
% axisX = 0:Re_x:Num_Cell_x;
% axisY = 0:Re_y:Num_Cell_y; %����Ϊ�ֱ���
% numX = length(axisX); %ȡ����
% numY = length(axisY);


%% ---------- initial distribution of target state 
delta_p=0.05;  % the stardard deviation of the inital distribution of position
delta_v=0.01;   % the stardard deviation of the inital distribution of velocity
new_velocity_Variance=delta_v^2;             % standard deviation 0.1;
% q2=2;                  % q2 the measurement noise
          
%% ---------- transition matrix          
F = [1 T_step 0   0; 
    0    1    0   0; 
    0    0    1   T_step; 
    0    0    0   1];            
           
Q=q1*[T_step^3/3  T_step^2/2  0           0 ;
        T_step^2/2  T_step      0           0 ;
        0             0         T_step^3/3  T_step^2/2;
        0             0         T_step^2/2  T_step];                         % ProcessNoise covariance matrix          
%Q = diag(20, 0.2, 20, 0.2);                      % JMPD��������ProcessNoise covariance matrix     
%% ---------- Rao-balckwellision parameters
Q_l=q1*[T_step,0;0,T_step];
Q_n=q1*[T_step^3/3,0;0,T_step^3/3];
Q_ln=q1*[T_step^2/2,0;0,T_step^2/2];
A_1_t=eye(2)-Q_ln/(Q_n)*T_step;
Q_1_l=Q_l-(Q_ln.')*(inv(Q_n))*Q_ln;

%% ------ load the continuious target trajectory 

Target_number=3;                        

velocity_init = 1; %
[initx,x] = GenerateTarget(Target_number,velocity_init,Num_Cell_x,Num_Cell_y,Total_time,F,Q);

x_dis = ceil(x(1,:,:)/Re_x)*Re_x; %�ֱܷ��Ŀ��λ�ã� ceil���������ȡ��
y_dis = ceil(x(3,:,:)/Re_y)*Re_y;

%% ---------- MC parameters            
repeati = 37;
SNR_T = [6,9,12];
SNR_num = length(SNR_T);  %length()������������ĳ���

% NpT = 2.^[6,8,9,10];
Np_T = 512;
Np_num = length(Np_T);
    
E_target_state_MC=zeros(7,Total_time,Target_number,repeati,SNR_num,Np_num); %��7��ʱ�䣬MC�����仯��SNR�仯���������仯��
single_run_time=zeros(Total_time,repeati,SNR_num,Np_num);
% figure(12);hold on;plot(squeeze(Pre_T_particle(1,t,:,:)),squeeze(Pre_T_particle(4,t,:,:)),'c.',x_c(1:3:13,t),x_c(2:3:14,t),'rx','LineWidth',2,'MarkerSize',8);

%% ===================== Particle filtering =========================
for Np_i=1:Np_num%Np_num%3%  %��ʾ����Np_num�β�ͬ�������������˲�
    Np = Np_T(Np_i);   %ÿ�ε�������       
    for SNR_i=1:SNR_num  %����SNR_num�β�ͬ����ȵ������˲�
        %����SNR����Ŀ����ȣ��õ�ÿһ֡�۲�ֵ
        SNR_dB=SNR_T(SNR_i);  %ÿ�ε������
        Frame_data = zeros(Num_Cell_y,Num_Cell_x,Total_time);
        Sigma_noise = 1;
        Frame_data = GenerateFrame(Num_Cell_y,Num_Cell_x,Total_time,SNR_dB,Sigma_noise,x,x_dis,y_dis,Target_number);
        for MC_i=1:repeati  %���ؿ��޴���
            display(['NP=',num2str(Np),'; SNR=',num2str(SNR_dB),'; MC=',num2str(MC_i)]);           


            %% ============== PF implementation ================   
            Pre_T_particle=zeros(7,Total_time,Target_number,Np);            % Pre_track particles
            Pre_T_particle_ori=zeros(7,Total_time,Target_number,Np);        % Pre_track particles
            Pre_T_life_index=zeros(Target_number,Np);                   % life length of the pre-tracks
            Pre_T_life_quality=zeros(Target_number,Np);                 % quality of life length of the pre-tracks
            Pre_weight0=zeros(Np,Total_time);               % Particles weights of Pre-track PF
            Pre_weight=zeros(Np,Total_time);                % Normolized Particles weights
            Pre_w_likehood_all=zeros(Np,Total_time);        % likehood part of the pre-PF weights
            Pre_track_bias=zeros(Np,Total_time);            % weight bias part of the pre-PF weights

   
            for t = 1:Time_point
                display(['Np=',num2str(Np_T(Np_i)),'; SNR=',num2str(SNR_T(SNR_i)),'; MC=',num2str(MC_i),'; t=',num2str(t)]);
                singerun_start=tic;                
                %% --------------- detection procedure ----------------    
                Detection_frame=Frame_data(:,:,t);    %�ѵ�ǰʱ�̵�����ƽ����Ϣ���浽����ʱ��ά�ȵĶ�λ���鷽�����
                Pre_track_Z=zeros(Target_number,Np);
                Partition_likehood=zeros(Target_number,Np);

                if t==1   
                    %disp(['t==1']);
                    Num_n_target=Target_number; % number of the birth pre-tracks     ����������Ŀ��ĸ���
                    index_x=initx(1,1:Target_number); %
                    index_y=initx(3,1:Target_number); %
                    index_vx=initx(2,1:Target_number); %
                    index_vy=initx(4,1:Target_number);   %initx������11��Ŀ��ĳ�ʼ״̬
                    % -------generate the new partitions of particles
                    %--------generate position based the detection measurements
                    position_x_p=repmat(index_x,Np,1)+delta_p*randn(Np,Num_n_target);
                    position_y_p=repmat(index_y,Np,1)+delta_p*randn(Np,Num_n_target);  
                    %--------generate velocity based on the detections
                    velocity_x_p=repmat((index_vx),Np,1);
                    velocity_y_p=repmat((index_vy),Np,1);  %��֪������Ϣ��������
                    %--------generate velocity variance
                    velocity_p_kk1=new_velocity_Variance.*ones(Np,2*Num_n_target);

                    %--------new_pretrack=zeros(4,Num_n_target,Np);
                    Pre_T_life_index=ones(Num_n_target,Np);
                    Pre_T_life_quality=ones(Num_n_target,Np);
                    for i=1:Np
                        Pre_T_particle(1:6,t,:,i)=[position_x_p(i,:);velocity_x_p(i,:);velocity_p_kk1(i,1:Num_n_target);position_y_p(i,:);velocity_y_p(i,:);velocity_p_kk1(i,1:Num_n_target)];
                    end   %�������ӵĸ�����Ϣ��λ�ã��ٶȣ����������ٶ�˥����
                    Pre_T_particle(7,t,:,:)=1;

                    particle_likehood_after_baise=ones(Target_number,Np);   %��ʼȨֵ��Ϊ1��
                else

                    %% --------------- evolution of the pre-tracks ----------------
                    %% -----------independent partition particle filter------------
                    Pre_track_Z=zeros(Target_number,Np);
                    Partition_likehood=zeros(Target_number,Np);
                    particle_likehood_after_baise=zeros(Target_number,Np);

                    for i=1:Target_number% 1%
                        for j=1:Np % 6%                  
                            %% === Rao-blackwellisation 
                            Pre_T_particle(1:6,t,i,j)= sample_RB( Pre_T_particle(1:6,t-1,i,j),T_step,Q_l,Q_n,Q_ln,A_1_t,Q_1_l,q1 ); %����������  Pre_T_particle(1:6,t-1,i,j)����һ��ʱ�̵�����
                            %%       %Pre_T_particle(1:6,t,i,j)��ʾtʱ�̵�i��Ŀ���j�����ӵ�6����Ϣ
                            Pre_T_life_index(i,j)=Pre_T_life_index(i,j)+1; 
                            Z_x_index=ceil(Pre_T_particle(1,t,i,j));
                            Z_y_index=ceil(Pre_T_particle(4,t,i,j));
                            if Z_x_index<=Num_Cell_x && Z_x_index>0 && Z_y_index<=Num_Cell_y && Z_y_index>0   %�ж��Ƿ���ƽ�淶Χ��                   
                                Pre_track_Z(i,j)=Detection_frame(Z_y_index,Z_x_index);
                                Pre_T_life_quality(i,j)=Pre_T_life_quality(i,j)+Detection_frame(Z_y_index,Z_x_index);
                                Pre_T_particle(7,t,i,j)=Detection_frame(Z_y_index,Z_x_index);
                                %% Gaussian likelihood ratio
                                %Partition_likehood(i,j)=exp(0.5*(2*Detection_frame(Z_y_index,Z_x_index)*Signal_amptitude-Signal_amptitude^2));
                                %% Rayleigh likelihood ratio or just likelihood
                                Partition_likehood(i,j)=raylpdf(Detection_frame(Z_y_index,Z_x_index),sqrt(Sigma_noise+Signal_amptitude^2))./raylpdf(Detection_frame(Z_y_index,Z_x_index),Sigma_noise);
                            else
                                Partition_likehood(i,j)=0;
                            end
                        end
                        Partition_likehood(i,:)=Partition_likehood(i,:)./sum(Partition_likehood(i,:));  %��һ��
                        %% === sample index funciton
                        [index_sample]=Sample_index(Partition_likehood(i,:));   %�ز���ѡ������ӵı��  ���ﲻ���ز�����
                        Pre_T_particle(:,:,i,:)=Pre_T_particle(:,:,i,index_sample);   %������ѡ����
                        Pre_T_life_quality(i,:)=Pre_T_life_quality(i,index_sample);   %�����Ӧ��Ȩֵ��
                        %% === retain the bias of sample: the likehood 
                        particle_likehood_after_baise(i,:)=Partition_likehood(i,index_sample);
                    end       
                end

                %% ---------- weights calculate of the Pre-track PF --------------- 
                cc=zeros(Np,1);
                for pre_Np_i=1:Np
                    %% ----sensor model: likelihood calculation                    
                    position_in=ceil(squeeze(Pre_T_particle([1,4],t,:,pre_Np_i)));
                    [ Pre_w_likehood_all(pre_Np_i,t)] = likelihood_calc( Detection_frame,position_in,Signal_amptitude);
                    %% -----Independent PF weights biase          
                    Pre_track_bias(pre_Np_i,t)=prod(particle_likehood_after_baise(:,pre_Np_i));
                    %% ------calculate the weights
                    Pre_weight0(pre_Np_i,t)= Pre_w_likehood_all(pre_Np_i,t)/Pre_track_bias(pre_Np_i,t);   
                end
                %% ------------- Normalize th Pre-weights -----------
                Pre_weight(:,t)=Pre_weight0(:,t)./sum(Pre_weight0(:,t));

                %% ------------ Resampling of the Pre-track PF ----------
                inIndex=1:Np;
                outIndex = deterministicR(inIndex,Pre_weight(:,t));
                Pre_T_particle_ori=Pre_T_particle;            % particles before resampling
                Pre_T_particle(:,:,:,:)=Pre_T_particle(:,:,:,outIndex);
                % Pre_T_life_index is not change during the resampling
                Pre_T_life_quality(:,:)=Pre_T_life_quality(:,outIndex); 
                
                %% ------------ time recording ----------
                single_run_time(t,MC_i,SNR_i,Np_i)=toc(singerun_start);%

%                 keyboard;
            end
            %% record the estimates
            E_target_state=zeros(7,Total_time,Target_number);
            for ii_t=1:Target_number
                E_target_state(:,:,ii_t)=mean(Pre_T_particle(:,:,ii_t,:),4);
            end
            E_target_state_MC(:,:,:,MC_i,SNR_i,Np_i)=E_target_state;
        end
    end      
end

%% =======================================================
%% ======== measurement calculation and plot =============
%% RMSE caculation not allowing swap
 