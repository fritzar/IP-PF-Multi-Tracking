function [ likelihood_out ] = likelihood_calc( Detection_frame,particle_in,Signal_amptitude)
%Pre_T_particle(:,t,i,pre_Np_i)
%LIKELIHOOD_CALC ¼ÆËãËÆÈ»º¯Êý
%   Detailed explanation goes here
format long;
Distinct_cell_n=0;            % the number of distinct cells  
Existing_cell=zeros(2,0);     % the position of each distinct cells
Existing_cell_tn=zeros(1,0);  % multiplicities
likelihood_out=0;
Target_number=size(particle_in,2);
Num_Cell_x=size(Detection_frame,2);
Num_Cell_y=size(Detection_frame,1);
Sigma_noise=1;
%% find out the number of distinct cells and the multiplicities of the
%% instinc cells
for i_t=1:Target_number
    %particle_temp=particle_in(:,i_t);
    x_temp=ceil(squeeze(particle_in(1,i_t)));                
    y_temp=ceil(squeeze(particle_in(2,i_t)));  
    
    if x_temp<=Num_Cell_x && x_temp>0 && y_temp<=Num_Cell_y && y_temp>0
        JD=sum(abs(Existing_cell-repmat([x_temp;y_temp],1,Distinct_cell_n)),1);
        if isempty(find(JD==0, 1)) || isempty(JD)
            Existing_cell = cat(2, Existing_cell, [x_temp;y_temp]);
            Distinct_cell_n=Distinct_cell_n+1;
            Existing_cell_tn(1,Distinct_cell_n)=1;
        else
            index_same = find(JD==0);
            Existing_cell_tn(1,index_same)=Existing_cell_tn(1,index_same)+1;
        end
    end  
end

%% likelihood calculation  %ÀÛ³Ë
if Distinct_cell_n>0
    likelihood_out=1;
    for i_dt=1:Distinct_cell_n
        x_temp=Existing_cell(1,i_dt);                
        y_temp=Existing_cell(2,i_dt); 
        %% Gassian likelihood ratio
%         likelihood_temp=exp(0.5*(2*Detection_frame(y_temp,x_temp)*Existing_cell_tn(1,i_dt)*Signal_amptitude-(Existing_cell_tn(1,i_dt)*Signal_amptitude)^2));
        %% Rayleigh likelihood or ratio                            
        likelihood_temp=raylpdf(Detection_frame(y_temp,x_temp),sqrt(Sigma_noise+Existing_cell_tn(1,i_dt)*Signal_amptitude^2))./raylpdf(Detection_frame(y_temp,x_temp),Sigma_noise);

        likelihood_out=likelihood_out*likelihood_temp;
    end
end

end

