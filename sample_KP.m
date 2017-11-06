function [ Particle_next ] = sample_KP( Particle,F,Q)
%SAMPLE_RB Summary of this function goes here
%   直接带入KP
Particle_next = zeros(1,6);
temp = [Particle(1:2) Particle(4:5)];
Particle_temp = reshape(temp,4,1);

processNoise = sample_gaussian(zeros(length(Q),1),Q,1)'; %each time calculate one frame
Particle_temp(1:4)=F*Particle_temp(1:4)+processNoise;  %选取动态转移概率作为重要性采样密度

Particle_temp = reshape(Particle_temp,1,4);

% Particle_next(1:2) = Particle_temp(1:2);
% Particle_next(4:5) = Particle_temp(3:4);
% 
% Particle_next(3) = Particle(3);
% Particle_next(6) = Particle(6);
Particle_next = [Particle_temp(1:2),Particle(3),Particle_temp(3:4), Particle(6)];
end