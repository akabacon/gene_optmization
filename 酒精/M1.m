function cost = M1(parameters)
    % 解析参数
    k = parameters(1);
    q = parameters(2);
    r = parameters(3);

    % 已知的时间与酒精含量数据
    time = [0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
    alcohol_content = [30, 68, 75, 82, 82, 77, 68, 68, 58, 51, 50, 41, 38, 35, 28, 25, 18, 15, 12, 10, 7, 7, 4];

    % 使用指定模型计算预测值
    predicted_alcohol_content = k * (exp(1).^(q .* time) - exp(1).^(r .* time));

    % 计算预测值与实际值之间的误差
    cost = sum((predicted_alcohol_content^2 - alcohol_content^2).);
    fprintf("k:%d q:%d r:%d\n",k,q,r);
end
