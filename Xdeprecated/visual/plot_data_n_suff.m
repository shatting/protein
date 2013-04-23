function [ ] = plot_data_n_suff( data )
%PLOT_DATA_N_SUFF Summary of this function goes here
%   Detailed explanation goes here

plot(data(:,1),data(:,2),'.'); 
suff_plot(suffstat_easy(data),[0,0,0],'mean');
