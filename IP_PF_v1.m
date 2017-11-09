% ==================
clear all;
close all;
clc;
format long; 

%% ======The initialization of the basic parameters=====
colorParticle={'bo-.','r+-.','k*-.','g>';'g^','k^','b^','y^';'bo-','ro','mo','go'};
Num_Cell_x=50;
Num_Cell_y=50;

Total_time=30;     % the number of total integrated frames
Total_time_data=30;
Time_point=30;
T_step=1;          % The size of the time cell:Time_step

q1=0.001;          % q1,q2 the level of process noise in target motion

Re_x = 1; %Resolution of x-axis
Re_y = 1;
% axisX = 0:Re_x:Num_Cell_x;
% axisY = 0:Re_y:Num_Cell_y; %步长为分辨率
% numX = length(axisX); %取长度
% numY = length(axisY);


%% ---------- initial distribution of target state 
delta_p=20;  % the stardard deviation of the inital distribution of position
delta_v=5;   % the stardard deviation of the inital distribution of velocity
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
%Q = diag(20, 0.2, 20, 0.2);                      % JMPD文中所用ProcessNoise covariance matrix     
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
[initx,x] = Generate_IP_Target(Target_number,velocity_init,Num_Cell_x,Num_Cell_y,Total_time,F,Q);

x_dis = ceil(x(1,:,:)/Re_x)*Re_x; %能分辨的目标位置， ceil朝正无穷方向取整
y_dis = ceil(x(3,:,:)/Re_y)*Re_y;

%% ---------- MC parameters            
repeati = 30;
SNR_T = [6];
SNR_num = length(SNR_T);  %length()函数，求数组的长度

% NpT = 2.^[6,8,9,10];
Np_T = 500;
Np_num = length(Np_T);
    
E_target_state_MC=zeros(7,Total_time,Target_number,repeati,SNR_num,Np_num); %（7，时间，MC次数变化，SNR变化，粒子数变化）
single_run_time=zeros(Total_time,repeati,SNR_num,Np_num);

Target_p_error = zeros(Total_time,repeati,Target_number,SNR_num,Np_num);
error_P = zeros(Total_time,Target_number,SNR_num,Np_num);

% figure(12);hold on;plot(squeeze(Pre_T_particle(1,t,:,:)),squeeze(Pre_T_particle(4,t,:,:)),'c.',x_c(1:3:13,t),x_c(2:3:14,t),'rx','LineWidth',2,'MarkerSize',8);

%% ===================== Particle filtering =========================
for Np_i=1:Np_num%Np_num%3%  %表示进行Np_num次不同粒子数的粒子滤波
    Np = Np_T(Np_i);   %每次的粒子数       
    for SNR_i=1:SNR_num  %进行SNR_num次不同信噪比的粒子滤波
        %根据SNR加上目标幅度，得到每一帧观测值
        SNR_dB=SNR_T(SNR_i);  %每次的信噪比
        Signal_amplitude=(10.^(SNR_dB./10)).^0.5;
        Frame_data = zeros(Num_Cell_y,Num_Cell_x,Total_time);
        Sigma_noise = 1;
        Frame_data = GenerateFrame(Num_Cell_y,Num_Cell_x,Total_time,SNR_dB,Sigma_noise,x,x_dis,y_dis,Target_number);
        for MC_i=1:repeati  %蒙特卡罗次数
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
                Detection_frame=Frame_data(:,:,t);    %把当前时刻的数据平面信息保存到不带时间维度的二位数组方便调用
                Pre_track_Z=zeros(Target_number,Np);
                Partition_likehood=zeros(Target_number,Np);

                if t==1   
                    %disp(['t==1']);
                    Num_n_target=Target_number; % number of the birth pre-tracks     航迹数就是目标的个数
                    index_x=initx(1,1:Target_number); %
                    index_y=initx(3,1:Target_number); %
                    index_vx=initx(2,1:Target_number); %
                    index_vy=initx(4,1:Target_number);   %initx保存了各个目标的初始状态
                    % -------generate the new partitions of particles
                    %--------generate position based the detection measurements
                    position_x_p=repmat(index_x,Np,1)+delta_p*randn(Np,Num_n_target);
                    position_y_p=repmat(index_y,Np,1)+delta_p*randn(Np,Num_n_target);
                    %% --------初始粒子均匀分布
%                     position_x_p = random('unif',1,50,Np,Target_number);
%                     position_y_p = random('unif',1,50,Np,Target_number);
                    %--------generate velocity based on the detections
                    velocity_x_p=repmat((index_vx),Np,1);
                    velocity_y_p=repmat((index_vy),Np,1);  %已知先验信息产生粒子
                    %--------generate velocity variance
                    velocity_p_kk1=new_velocity_Variance.*ones(Np,2*Num_n_target);

                    %--------new_pretrack=zeros(4,Num_n_target,Np);
                    Pre_T_life_index=ones(Num_n_target,Np);
                    Pre_T_life_quality=ones(Num_n_target,Np);
                    for i=1:Np
                        Pre_T_particle(1:6,t,:,i)=[position_x_p(i,:);velocity_x_p(i,:);velocity_p_kk1(i,1:Num_n_target);position_y_p(i,:);velocity_y_p(i,:);velocity_p_kk1(i,1:Num_n_target)];
                    end   %保存粒子的各种信息，位置，速度，速度波动
                    Pre_T_particle(7,t,:,:)=1;

                    particle_likehood_bias=ones(Target_number,Np);   %初始权值都为1
                else

                    %% --------------- evolution of the pre-tracks ----------------
                    %% -----------independent partition particle filter------------
                    Pre_track_Z=zeros(Target_number,Np);
                    Partition_likehood=zeros(Target_number,Np);
                    particle_likehood_bias=zeros(Target_number,Np);

                    for i=1:Target_number% 1）For each partition
                        for j=1:Np % 6%                  
                            %% === Rao-blackwellisation 
                            %Pre_T_particle(1:6,t,i,j)= sample_RB( Pre_T_particle(1:6,t-1,i,j),T_step,Q_l,Q_n,Q_ln,A_1_t,Q_1_l,q1 ); %产生新粒子  Pre_T_particle(1:6,t-1,i,j)是上一个时刻的粒子
                            Pre_T_particle(1:6,t,i,j)= sample_KP(Pre_T_particle(1:6,t-1,i,j),F,Q); %纯粹的粒子滤波
                            %%       %Pre_T_particle(1:6,t,i,j)表示t时刻第i个目标第j个粒子的6个信息
                            Pre_T_life_index(i,j)=Pre_T_life_index(i,j)+1; 
                            Z_x_index=ceil(Pre_T_particle(1,t,i,j));
                            Z_y_index=ceil(Pre_T_particle(4,t,i,j));
                            if Z_x_index<=Num_Cell_x && Z_x_index>0 && Z_y_index<=Num_Cell_y && Z_y_index>0   %判断是否在平面范围内                   
                                Pre_track_Z(i,j)=Detection_frame(Z_y_index,Z_x_index);
                                Pre_T_life_quality(i,j)=Pre_T_life_quality(i,j)+Detection_frame(Z_y_index,Z_x_index);
                                Pre_T_particle(7,t,i,j)=Detection_frame(Z_y_index,Z_x_index);
                                %% Gaussian likelihood ratio
                                %Partition_likehood(i,j)=exp(0.5*(2*Detection_frame(Z_y_index,Z_x_index)*Signal_amplitude-Signal_amplitude^2));
                                %% Rayleigh likelihood ratio or just likelihood
                                Partition_likehood(i,j)=raylpdf(Detection_frame(Z_y_index,Z_x_index),sqrt(Sigma_noise+Signal_amplitude^2))./raylpdf(Detection_frame(Z_y_index,Z_x_index),Sigma_noise);
                            else
                                Partition_likehood(i,j)=0;
                            end
                        end
                        Partition_likehood(i,:)=Partition_likehood(i,:)./sum(Partition_likehood(i,:));  %归一化
                        %% === sample index funciton
                        [index_sample]=Sample_index(Partition_likehood(i,:));   %重采样选择的粒子的编号  
                        Pre_T_particle(:,:,i,:)=Pre_T_particle(:,:,i,index_sample);   %复制所选粒子
                        Pre_T_life_quality(i,:)=Pre_T_life_quality(i,index_sample);   %保存对应的权值
                        %% === retain the bias of sample: 保存观测信息用于优化重要性密度函数 
                        particle_likehood_bias(i,:)=Partition_likehood(i,index_sample); 
                    end       
                end

                %% ---------- weights calculate of the Pre-track PF --------------- 
                cc=zeros(Np,1);
                for pre_Np_i=1:Np
                    %% ----sensor model: likelihood calculation                    
                    position_in=ceil(squeeze(Pre_T_particle([1,4],t,:,pre_Np_i)));
                    [ Pre_w_likehood_all(pre_Np_i,t)] = likelihood_calc( Detection_frame,position_in,Signal_amplitude);
                    %% -----Independent PF weights biase          
                    Pre_track_bias(pre_Np_i,t)=prod(particle_likehood_bias(:,pre_Np_i));
                    %% ------calculate the weights
                    Pre_weight0(pre_Np_i,t)= Pre_w_likehood_all(pre_Np_i,t)/Pre_track_bias(pre_Np_i,t);   
                end
                %% ------------- Normalize the Pre-weights -----------
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

%                keyboard;
            end
            %% record the estimates
            E_target_state=zeros(7,Total_time,Target_number);
            for t_i=1:Target_number
                E_target_state(:,:,t_i)=mean(Pre_T_particle(:,:,t_i,:),4);
            end
            E_target_state_MC(:,:,:,MC_i,SNR_i,Np_i)=E_target_state;
        end
        
        for t_i= 1:Target_number
            Target_p_error(:,:,t_i,SNR_i,Np_i) = (squeeze(E_target_state_MC(1,:,t_i,:,SNR_i,Np_i))-repmat(x(1,((t_i-1)*Total_time+1):(t_i*Total_time)),repeati,1)').^2 + (squeeze(E_target_state_MC(4,:,t_i,:,SNR_i,Np_i))-repmat(x(1,((t_i-1)*Total_time+1):(t_i*Total_time)),repeati,1)').^2; %注意对真实目标的矩阵有转置！T*Monte
            error_P(:,t_i,SNR_i,Np_i) = sqrt(mean(Target_p_error(:,:,t_i,SNR_i,Np_i),2)); %%T，各帧的RMSE
        end
    end      
end

%% =======================================================
%% ======== measurement calculation and plot =============
figure(2)
hold on
for n = 1:Target_number
    plot(squeeze(x(1,n,:)),squeeze(x(3,n,:)),colorParticle{1,n},'Linewidth',2)
    plot(E_target_state(1,:,n),E_target_state(4,:,n),'g^:','Linewidth',3)
end
axis([0,50,0,50])
grid on
legend ('Target1','Estimation1','Target2','Estimation2','Target3','Estimation3')
title('三目标跟踪结果')
xlabel('x方向距离')
ylabel('y方向距离')
% 
%% RMSE caculation not allowing swap
figure(60)
hold on;
plot(error_P(:,1,1,1),'b^-');
plot(error_P(:,2,1,1),'r^-');
plot(error_P(:,3,1,1),'k^-');
title('三个目标的各帧均方误差')
% axis([0,Total_time,0,1])
xlabel('时间/帧')
ylabel('均方误差')
legend('Target1','Target2','Target3')
grid on
hold off

RMS_t = zeros(Target_number);
RMS_t(:) = mean(error_P(:,1,1,1),1)


 