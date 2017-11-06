function [ index_sample ] = Sample_index( Partition_likehood)
%SAMPLE_INDEX Summary of this function goes here
%   Detailed explanation goes here
% Partition_likehood=weight_1st;
index_sample=zeros(1,length(Partition_likehood));
CDF_weights=zeros(1,length(Partition_likehood));
CDF_weights(1)=Partition_likehood(1);
for i=2:length(Partition_likehood)
    CDF_weights(i)=CDF_weights(i-1)+Partition_likehood(i);   %CDF_weights(i)保存的是从1到i的似然的累加
end
% CDF_weights
B=([1:length(Partition_likehood)]-rand(1))./length(Partition_likehood);

j=1;
for i=1:length(Partition_likehood)
    while B(i)>CDF_weights(j)
        j=j+1;
    end
    index_sample(i)=j;
end 

