clear; clc; close all;
VA = 15;    % A船航速 kn
VB = 14;    % B船航速 kn
S0 = 10;    % 初始距离 n mile
D0 = 1;     % 安全CPA阈值1 n mile，B船圆形半径
maxIter = 150;  % 迭代次数150次
popSize = 50;   % 种群总数
selectNum = 25; % 轮盘赌选出25个父代
crossNum = 25;  % 交叉变异生成25个子代

% 变量上下限 [t1(min), t2(min), a1(°), a2(°)]
lb = [0,    3,   5,   8];
ub = [18,  20,  45,  75];

fitnessFun = @(x) calcFitness(x,VA,VB,S0,D0);

pop = rand(popSize,4).*(ub-lb) + lb;
fitRecord = zeros(maxIter,1);

for iter = 1:maxIter
    fitVal = zeros(popSize,1);
    for i = 1:popSize
        fitVal(i) = fitnessFun(pop(i,:));
    end
    
    [bestFit,bestIdx] = max(fitVal);
    fitRecord(iter) = bestFit;
    bestInd = pop(bestIdx,:);
    
    fitShift = fitVal - min(fitVal) + 1e-6;
    prob = fitShift / sum(fitShift);
    parentIdx = zeros(selectNum,1);
    for s = 1:selectNum
        r = rand();
        cumsumProb = 0;
        for i = 1:popSize
            cumsumProb = cumsumProb + prob(i);
            if r <= cumsumProb
                parentIdx(s) = i;
                break;
            end
        end
    end
    parentPop = pop(parentIdx,:);
    
    childPop = zeros(crossNum,4);
    for i = 1:crossNum
        p1 = parentPop(randi(selectNum),:);
        p2 = parentPop(randi(selectNum),:);
        crossPos = randi(3);
        child = [p1(1:crossPos), p2(crossPos+1:end)];
        if rand < 0.1
            child = child + randn(1,4).*0.5;
        end
        child = max(lb, min(ub, child));
        childPop(i,:) = child;
    end
    pop = [parentPop; childPop];
end

fprintf('==================== 最优避让参数 ====================\n');
t1_opt = bestInd(1);
t2_opt = bestInd(2);
a1_opt = bestInd(3);
a2_opt = bestInd(4);
fprintf('直航时间 t1 = %.2f min\n',t1_opt);
fprintf('避让航行 t2 = %.2f min\n',t2_opt);
fprintf('避让转角 a1 = %.2f °\n',a1_opt);
fprintf('复航转角 a2 = %.2f °\n',a2_opt);

[CPA_opt,TCPA_opt] = calcCPA_TCPA(bestInd,VA,VB,S0);
fprintf('CPA = %.4f n mile (安全阈值D0=1)\n',CPA_opt);
fprintf('TCPA = %.2f min\n',TCPA_opt);

figure('Name','GA收敛曲线');
iterAxis = 1:maxIter;
winSize = 8;
fitSmooth = movmean(fitRecord, winSize);
iterFine = linspace(1, maxIter, 800);
fitFine = spline(iterAxis, fitSmooth, iterFine);
plot(iterFine, fitFine, 'b-','LineWidth',2.2);
xlabel('迭代次数','FontSize',11); 
ylabel('当代最优适应度','FontSize',11);
title('遗传算法适应度收敛曲线','FontSize',13);
grid on; grid minor;
legend('平滑收敛曲线','Location','best');

dynamicShipAnim(bestInd,VA,VB,S0,D0);

function [CPA,TCPA] = calcCPA_TCPA(x,VA,VB,S0)
    t1 = x(1); t2 = x(2); a1 = x(3);
    rad_a1 = deg2rad(a1);
    CPA = (t2 * sin(rad_a1)) / 4;
    S1 = S0 - (VA+VB)*t1/60;
    Vrx = VA*cos(rad_a1) + VB;
    dt = (60 * S1) / Vrx;
    TCPA = t1 + dt;
end

function fit = calcFitness(x,VA,VB,S0,D0)
    t1 = x(1); t2 = x(2); a1 = x(3); a2 = x(4);
    rad1 = deg2rad(a1); rad2 = deg2rad(a2);
    totalT = t1 + t2 + 12;
    dtStep = 0.2;
    tList = 0:dtStep:totalT;
    backBaseTime = inf;
    
    [CPA,TCPA] = calcCPA_TCPA(x,VA,VB,S0);
    penaltyCPA = 0;
    if CPA < D0
        penaltyCPA = 200*(D0 - CPA);
    end
    
    % 全程任意时刻船距小于安全阈值
    penaltyDist = 0;
    for t = tList
        [Ax,Ay] = calcShipApos(t,t1,t2,a1,a2,VA);
        Bx = S0 - VB * t / 60;
        By = 0;
        dist_now = sqrt((Ax-Bx)^2 + (Ay-By)^2);
        if dist_now < D0
            penaltyDist = penaltyDist + 300*(D0 - dist_now);
        end
        % 记录第一次回到基线y≈0的时刻
        if abs(Ay) < 0.01 && t > t1 + t2
            backBaseTime = min(backBaseTime, t);
        end
    end
    
    penaltyBackBase = 0;
    if backBaseTime == inf
        penaltyBackBase = 250;
    end
    
    rad_a1 = deg2rad(a1);
    S1 = S0 - (VA+VB)*t1/60;
    Vrx = VA*cos(rad_a1)+VB;
    dt = 60*S1 / Vrx;
    penaltyTime = 0;
    if dt > t2 || S1 < 0
        penaltyTime = 15;
    end
    penaltyAngle = 0;
    if a2 <= a1
        penaltyAngle = 10;
    end
    
    objCPA = -abs(CPA - D0);
    objTCPA = -abs(TCPA - 15);
    objT1 = t1 / 18;
    totalObj = 0.8*objCPA + 0.15*objTCPA + 0.05*objT1;
    
    fit = totalObj - penaltyCPA - penaltyDist - penaltyBackBase - penaltyTime - penaltyAngle;
end

function [Ax,Ay] = calcShipApos(t,t1,t2,a1,a2,VA)
    rad1 = deg2rad(a1); rad2 = deg2rad(a2);
    persistent backFlag; % 持久标记：是否已经回到基线
    if t == 0
        backFlag = 0; % 仿真初始化重置标记
    end
    
    if t <= t1
        % 阶段1：直航基线
        Ax = VA * t / 60;
        Ay = 0;
        backFlag = 1; % 初始就在基线上
    elseif t <= t1 + t2
        % 阶段2：避让转向
        dtSeg = t - t1;
        x0 = VA * t1 / 60;
        Ax = x0 + VA * dtSeg /60 * cos(rad1);
        Ay = VA * dtSeg /60 * sin(rad1);
        backFlag = 0; % 离开基线，标记清零
    else
        dt3 = t - t1 - t2;
        x0 = VA*t1/60 + VA*t2/60*cos(rad1);
        y0 = VA*t2/60*sin(rad1);
        Ax_temp = x0 + VA * dt3 /60 * cos(rad1 - rad2);
        Ay_temp = y0 + VA * dt3 /60 * sin(rad1 - rad2);
        
        % 只要标记已触发，永久锁定基线
        if backFlag == 1
            Ax = x0 + VA * dt3 / 60;
            Ay = 0;
        else
            % 未回基线：正常复航，同时检测是否触达基线
            if abs(Ay_temp) < 0.01
                backFlag = 1; % 首次触达，永久锁定
                Ax = x0 + VA * dt3 / 60;
                Ay = 0;
            else
                Ax = Ax_temp;
                Ay = Ay_temp;
            end
        end
    end
end

function dynamicShipAnim(x,VA,VB,S0,D0)
    t1 = x(1); t2 = x(2); a1 = x(3); a2 = x(4);
    dtStep = 0.08;
    totalT = t1 + t2 + 12;
    tList = 0:dtStep:totalT;
    
    figure('Name','双船避碰动态轨迹分析图');
    hold on; grid on; axis equal;
    xlabel('X 坐标 (n mile)');
    ylabel('Y 坐标 (n mile)');
    title('自动船对遇避让全过程动态仿真');
    
    Ax_record = []; Ay_record = [];
    Bx_record = []; By_record = [];
    
    shipA = scatter(0,0,30,'r','filled','DisplayName','A让路船');
    theta = linspace(0,2*pi,100);
    circleX = D0 * cos(theta);
    circleY = D0 * sin(theta);
    shipB_circle = plot(S0 + circleX, 0 + circleY, 'b-','LineWidth',1.2,'DisplayName',['B直航船安全圈(D0=',num2str(D0),'nmile)']);
    shipB_center = scatter(S0,0,20,'b','filled');
    
    trajA_line = plot(NaN,NaN,'r-','LineWidth',1.2);
    trajB_line = plot(NaN,NaN,'b-','LineWidth',1.2);
    
    infoText = text(0.02,0.95,'','Units','normalized','FontSize',10,...
        'BackgroundColor',[0.98,0.98,0.98]);
    legend('Location','best');
    
    for idx = 1:length(tList)
        t = tList(idx);
        % 调用坐标函数，自动处理归基线后直行逻辑
        [Ax,Ay] = calcShipApos(t,t1,t2,a1,a2,VA);
        Bx = S0 - VB * t / 60;
        By = 0;
        
        Ax_record = [Ax_record, Ax];
        Ay_record = [Ay_record, Ay];
        Bx_record = [Bx_record, Bx];
        By_record = [By_record, By];
        
        shipA.XData = Ax; shipA.YData = Ay;
        shipB_center.XData = Bx; shipB_center.YData = By;
        shipB_circle.XData = Bx + circleX;
        shipB_circle.YData = By + circleY;
        trajA_line.XData = Ax_record; trajA_line.YData = Ay_record;
        trajB_line.XData = Bx_record; trajB_line.YData = By_record;
        
        [CPA_now,TCPA_now] = calcCPA_TCPA(x,VA,VB,S0);
        dist_now = sqrt((Ax-Bx)^2 + (Ay-By)^2);
        strInfo = sprintf(['当前航行时间：%.2f min\n',...
            '安全CPA阈值：%.2f n mile\n',...
            '两船实时距离：%.3f n mile\n',...
            '本次避让CPA：%.3f n mile\n',...
            'TCPA：%.2f min'],t,D0,dist_now,CPA_now,TCPA_now);
        infoText.String = strInfo;
        
        drawnow;
        pause(0.01);
    end
    text(0.35,0.05,'避让仿真动画结束','FontSize',14,'Color','magenta','Units','normalized');
end