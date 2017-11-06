rmse=zeros(1,6);
rmse=[1.8998,1.2519,1.2028,1.0814,0.3704,0.4694];
SNR_t=[3,4,6,8,10,12];

%{
for i=1:6
    hold on;
    grid on;
    plot(SNR_t(i),rmse(i),'+');
    axis([0 13 0 2]);
    xlabel('SNR/dB');ylabel('RMSE');
end
%}
figure(10);
plot(SNR_t,rmse,'+r-');