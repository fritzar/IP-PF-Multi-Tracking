function [initx,x] = GenerateTarget(Target_number,velocity_init,Num_Cell_x,Num_Cell_y,Total_time_data,F,Q)
%生成目标
colorParticle={'bo-','r+:','g*-.','k.';'g^','k^','b^','y^';'bo-','ro','mo','go'};
initx=zeros(4,Target_number); %x=[x,x',y,y']T
x = zeros(4,Target_number,Total_time_data); %每Total_time_data个列存放一个目标的各帧状态
x_c = zeros(3*Target_number,Total_time_data); %why 3?
% F = [1 T_step 0 0
%      0   1    0 0
%      0   0    1 T_step
%      0   0    0 1];
% Q=[T_step^3*q1/3  T_step^2*q1/2  0           0 ;
%    T_step^2*q1/2 T_step*q1      0           0 ;
%    0             0           T_step^3*q1/3  T_step^2*q1/2;
%    0             0           T_step^2*q1/2  q1*T_step];        % ProcessNoise covariance matrix
Sigma_noise=1; %sigma_noise
% 
% i_p = [10,6;6,Num_Cell_y/2;Num_Cell_x/2,Num_Cell_y-6;Num_Cell_x-6,10]'; %position i_p(4*2)? 
% v1=(0.8+0.4*rand(1,4))*velocity_init; 
% v = repmat(v1,2,1); %replicate v1 by 2*1
% alpha = [pi/4+0.2,0.1,-pi/2-0.1,3*pi/4+0.3];
% i_v = v.*[cos(alpha(1)),sin(alpha(1));cos(alpha(2)),sin(alpha(2));cos(alpha(3)),sin(alpha(3));cos(alpha(4)),sin(alpha(4))]'; %velocity
% initx([1,3],:) = i_p(:,1:Target_number);
% initx([2,4],:) = i_v(:,1:Target_number);

initx([1,3],1)  = [10;10];
initx([1,3],2) = [42.5,25];
initx([1,3],3) = [15;30];

initx([2,4],:) = [1,1;-1,0;1,-1]';



processNoise =Sigma_noise*sample_gaussian(zeros(length(Q),1),Q,Total_time_data)';% Note ' here! state_noise_samples is: ss-Total_time_allfor t=1:Total_time_all
for t = 1:Total_time_data
    for n = 1:Target_number
        if t==1
            x(:,n,t) = initx(:,n) + processNoise(t);
        else
            x(:,n,t) = F*x(:,n,t-1) + processNoise(t);
        end
    end
end

figure(1)
hold on
for n = 1:Target_number
    plot(squeeze(x(1,n,:)),squeeze(x(3,n,:)),colorParticle{1,n},'Linewidth',2)
end
hold off
axis([0,50,0,50])
title('三目标轨迹')
grid on
legend ('Target1','Target2(moving from right to left)','Target3')

