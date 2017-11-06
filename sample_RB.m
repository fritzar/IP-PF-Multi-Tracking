function [ Particle_T ] = sample_RB( Particle_O,T_step,Q_l,Q_n,Q_ln,A_1_t,Q_1_l,q1)
%SAMPLE_RB Summary of this function goes here
%   Detailed explanation goes here
        mean_temp=[Particle_O(1)+T_step*Particle_O(2);Particle_O(4)+T_step*Particle_O(5)];
        covar_temp=(T_step.^2).*[Particle_O(3),0;0, Particle_O(6)]+ [T_step^3*q1/3,0; 0,T_step^3*q1/3;];
        particle_temp=sample_gaussian(mean_temp,covar_temp,1);
        Particle_T(1)=particle_temp(1);
        Particle_T(4)=particle_temp(2);
        %end
     
     %% -- Kalman filter time update: prediction
            z_t=[Particle_T(1)-Particle_O(1);Particle_T(4)-Particle_O(4)];
            X_tt=[Particle_O(2);Particle_O(5)];
            P_tt=[Particle_O(3),0;0,Particle_O(6)];
            
            N_t=(T_step^2).*[Particle_O(3),0;0, Particle_O(6)]+Q_n;
            L_t=A_1_t*P_tt*T_step*inv(N_t);            
            
            V_m_t1_t=A_1_t*X_tt+Q_ln.'*(inv(Q_n))*z_t+L_t*(z_t-T_step*X_tt);
            V_p_t1_t=A_1_t*P_tt*(A_1_t.')+Q_1_l-L_t*N_t*(L_t.');
            
            Particle_T([2,5])=V_m_t1_t;
            Particle_T(3)=V_p_t1_t(1,1);
            Particle_T(6)=V_p_t1_t(2,2);


end

