clear
clc
data = xlsread('新增表.xlsx');
N = length(data);
X = data(1:(N-1),1:3);
Y = data(2:N,4);
X_train = X(1:40,:);
Y_train = Y(1:40);
% x2= data(1,1:(length(data)-10));
%X=[1 13 1.5; 1.4 19 3; 1.8 25 1; 2.2 10 2.5;2.6 16 0.5; 3 22 2; 3.4 28 3.5; 3.5 30 3.7];
%Y=[0.330; 0.336; 0.294; 0.476; 0.209; 0.451; 0.482; 0.5];
choose=2;
fit_nonlinear_data(X_train, Y_train, X, Y, choose)