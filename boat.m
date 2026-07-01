%% 双船对遇避碰 - 遗传算法优化避让参数 t1,t2,a1,a2
clear; clc; close all;
%% ===================== 1. 固定全局参数 =====================
VA = 15;    % A船航速 kn
VB = 14;    % B船航速 kn
S0 = 10;    % 初始距离 n mile
D0 = 2;     % 安全CPA阈值 n mile
maxIter = 100;  % GA最大迭代次数
popSize = 50;   % 种群总数
selectNum = 25; % 轮盘选择父代数量
crossNum = 25;  % 交叉变异产生子代数量

% 变量上下限 [t1, t2, a1, a2]
% t1(min), t2(min), a1(°), a2(°)
lb = [0,    3,   5,   8];
ub = [18,  20,  45,  75];

%% ===================== 2. 适应度函数 =====================
fitnessFun = @(x) calcFitness(x,VA,VB,S0,D0);

%% ===================== 3. 遗传算法主循环 =====================
% 初始化种群
pop = rand(popSize,4).*(ub-lb) + lb;
fitRecord = zeros(maxIter,1); % 记录每代最优适应度

for iter = 1:maxIter
    % 计算所有个体适应度
    fitVal = zeros(popSize,1);
    for i = 1:popSize
        fitVal(i) = fitnessFun(pop(i,:));
    end
    
    % 记录当代最优
    [bestFit,bestIdx] = max(fitVal);
    fitRecord(iter) = bestFit;
    bestInd = pop(bestIdx,:);
    
    %% 轮盘赌选择 25个父代
    fitNorm = fitVal - min(fitVal) + 1e-6; % 避免0
    prob = fitNorm / sum(fitNorm);
    parentIdx = randsample(1:popSize, selectNum, true, prob);
    parentPop = pop(parentIdx,:);
    
    %% 交叉+变异生成25个子代
    childPop = zeros(crossNum,4);
    for i = 1:crossNum
        p1 = parentPop(randi(selectNum),:);
        p2 = parentPop(randi(selectNum),:);
        % 单点交叉
        crossPos = randi(3);
        child = [p1(1:crossPos), p2(crossPos+1:end)];
        % 高斯变异
        if rand < 0.1
            child = child + randn(1,4).*0.5;
        end
        % 限制在边界内
        child = max(lb, min(ub, child));
        childPop(i,:) = child;
    end
    
    % 合并父代+子代 形成新种群50个
    pop = [parentPop; childPop];
end

%% ===================== 4. 输出最优参数 =====================
fprintf('==================== 最优避让参数 ====================\n');
t1_opt = bestInd(1);
t2_opt = bestInd(2);
a1_opt = bestInd(3);
a2_opt = bestInd(4);
fprintf('直航时间 t1 = %.2f min\n',t1_opt);
fprintf('避让航行 t2 = %.2f min\n',t2_opt);
fprintf('避让转角 a1 = %.2f °\n',a1_opt);
fprintf('复航转角 a2 = %.2f °\n',a2_opt);

% 计算最优参数下CPA、TCPA
[CPA_opt,TCPA_opt] = calcCPA_TCPA(bestInd,VA,VB,S0);
fprintf('CPA = %.4f n mile (安全阈值2)\n',CPA_opt);
fprintf('TCPA = %.2f min\n',TCPA_opt);

%% ===================== 5. 绘制适应度收敛曲线 =====================
figure('Name','GA收敛曲线');
plot(1:maxIter, fitRecord,'LineWidth',1.5);
xlabel('迭代次数'); ylabel('当代最优适应度');
title('遗传算法适应度收敛曲线'); grid on;

%% ===================== 6. 绘制A、B两船完整航行轨迹 =====================
simShipTraj(bestInd,VA,VB,S0);

%% ===================== 子函数1：计算CPA TCPA =====================
function [CPA,TCPA] = calcCPA_TCPA(x,VA,VB,S0)
    t1 = x(1); t2 = x(2); a1 = x(3);
    rad_a1 = deg2rad(a1);
    % CPA
    CPA = (t2 * sin(rad_a1)) / 4;
    % TCPA
    S1 = S0 - (VA+VB)*t1/60;
    Vrx = VA*cos(rad_a1) + VB;
    dt = (60 * S1) / Vrx;
    TCPA = t1 + dt;
end

%% ===================== 子函数2：适应度计算 =====================
function fit = calcFitness(x,VA,VB,S0,D0)
    t1 = x(1); t2 = x(2); a1 = x(3); a2 = x(4);
    [CPA,TCPA] = calcCPA_TCPA(x,VA,VB,S0);
    
    % 约束1：CPA必须≥安全距离D0
    penaltyCPA = 0;
    if CPA < D0
        penaltyCPA = 100*(D0 - CPA);
    end
    % 约束2：最近会遇必须在避让阶段内 dt <= t2
    rad_a1 = deg2rad(a1);
    S1 = S0 - (VA+VB)*t1/60;
    Vrx = VA*cos(rad_a1)+VB;
    dt = 60*S1 / Vrx;
    penaltyTime = 0;
    if dt > t2 || S1 < 0
        penaltyTime = 80;
    end
    % 约束3：复航角度大于避让角度，保证能回航
    penaltyAngle = 0;
    if a2 <= a1
        penaltyAngle = 50;
    end
    
    % 多目标加权收益：
    % 1. CPA略大于2最好，过大扣分；2. TCPA适中，不要过小；3. t1尽量大（晚转向）
    objCPA = -abs(CPA - D0);
    objTCPA = -abs(TCPA - 15); % 期望TCPA接近15min
    objT1 = t1 / 18;           % 希望t1尽可能大
    
    totalObj = 0.6*objCPA + 0.3*objTCPA + 0.1*objT1;
    fit = totalObj - penaltyCPA - penaltyTime - penaltyAngle;
end

%% ===================== 子函数3：船舶轨迹仿真绘图 =====================
function simShipTraj(x,VA,VB,S0)
    t1 = x(1); t2 = x(2); a1 = x(3); a2 = x(4);
    rad1 = deg2rad(a1); rad2 = deg2rad(a2);
    dtStep = 0.1; % 仿真步长 min
    totalT = t1 + t2 + 15; % 额外多仿真15min看复航远离
    tList = 0:dtStep:totalT;
    
    % B船坐标：全程沿-x匀速直行
    Bx = zeros(size(tList)); By = zeros(size(tList));
    % 初始位置：A(0,0), B(S0, 0) 对遇
    for i = 1:length(tList)
        t = tList(i);
        Bx(i) = S0 - VB * t / 60;
        By(i) = 0;
    end
    
    % A船三段运动坐标
    Ax = zeros(size(tList)); Ay = zeros(size(tList));
    for i = 1:length(tList)
        t = tList(i);
        if t <= t1
            % 阶段1：直航向B，沿+x
            Ax(i) = VA * t / 60;
            Ay(i) = 0;
        elseif t <= t1 + t2
            % 阶段2：避让航向，偏右a1
            dtSeg = t - t1;
            x0 = VA * t1 / 60;
            Ax(i) = x0 + VA * dtSeg /60 * cos(rad1);
            Ay(i) = 0 + VA * dtSeg /60 * sin(rad1);
        else
            % 阶段3：复航偏转a2，逐步拉回y=0原航线
            dt1 = t1;
            dt2 = t2;
            dt3 = t - t1 - t2;
            x0 = VA*dt1/60 + VA*dt2/60*cos(rad1);
            y0 = VA*dt2/60*sin(rad1);
            Ax(i) = x0 + VA * dt3 /60 * cos(rad1 - rad2);
            Ay(i) = y0 + VA * dt3 /60 * sin(rad1 - rad2);
        end
    end
    
    % 绘图
    figure('Name','A/B船舶避碰轨迹');
    plot(Ax,Ay,'r-','LineWidth',1.5,'DisplayName','A避让船'); hold on;
    plot(Bx,By,'b-','LineWidth',1.5,'DisplayName','B直航船');
    scatter(Ax(1),Ay(1),40,'r','filled'); text(Ax(1)+0.2,Ay(1),'A起点');
    scatter(Bx(1),By(1),40,'b','filled'); text(Bx(1)-0.6,By(1),'B起点');
    xlabel('X (n mile)'); ylabel('Y (n mile)');
    title('最优参数下两船航行避碰轨迹');
    legend; grid on; axis equal;
end