function [beta, r]=fit_nonlinear_data(X, Y, X_true,Y_true, choose)
% Input: X 自变量数据(N, D)， Y 因变量(N, 1)，choose 1-regress, 2-nlinfit 3-lsqcurvefit
if choose==1
    X1=[ones(length(X(:, 1)), 1), X];
    [beta, bint, r, rint, states]=regress(Y, X1);
    % 多元线性回归
    % y=beta(1)+beta(2)*x1+beta(3)*x2+beta(4)*x3+...
    % beta—系数估计
    % bint—系数估计的上下置信界
    % r—残差
    % rint—诊断异常值的区间
    % states—模型统计信息
    Y_pred = beta(1) + beta(2)*X_true(:,1) + beta(3)*X_true(:,2) + beta(4)*X_true(:,3);
    % Y_pred = beta(1) + beta(2)*X_true(:,1);
    plot(Y_true,'linewidth',4)
    hold on
    plot(Y_pred,'linewidth',4)
elseif choose==2
    beta0=ones(7, 1);
    [beta, r, J]=nlinfit(X, Y, @myfun, beta0)
    % 非线性回归
    % beta—系数估计
    % r—残差
    % J—雅可比矩阵
    [Y_pred,delta]=nlpredci(@myfun, X_true, beta, r, 'Jacobian', J)
    % 非线性回归预测置信区间
    % Ypred—预测响应
    % delta—置信区间半角
    %plot(X(:, 1), Y, 'k.', X(:, 1), Ypred, 'r');
    plot(Y_true,'linewidth',4)
    hold on
    plot(Y_pred,'linewidth',4)
    %saveas(gcf,sprintf('非线性曲线拟合_1.jpg'),'bmp');
end
end
 
 
