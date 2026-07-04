%% ====================== 主程序：船舶避碰帕累托双目标遗传避让优化（动态动画版） ======================
clear; clc; close all;

%% 1. 交互式输入船舶基础参数
disp('========== 船舶避碰综合决策系统(COLREGs+帕累托遗传避让优化) ==========');
own_crs  = input('请输入本船航向(°，正北0°顺时针)：');
own_spd  = input('请输入本船航速(节 kn)：');
azi0     = input('请输入目标船相对方位角(雷达观测本船看目标，°)：');
R0       = input('请输入两船当前相对距离(海里 n mile)：');
tar_crs  = input('请输入目标船航向(°)：');
tar_spd  = input('请输入目标船航速(节 kn)：');
Dsafe    = input('请设置安全会遇距离阈值(开阔水域推荐2海里)：');

%% 2. 角度转换 & 速度分解 & 原始DCPA/TCPA计算
deg2rad = pi / 180;
th0 = own_crs * deg2rad;
tht = tar_crs * deg2rad;
az  = azi0 * deg2rad;

u0 = own_spd * sin(th0); v0 = own_spd * cos(th0);
ut = tar_spd * sin(tht); vt = tar_spd * cos(tht);
du = ut - u0; dv = vt - v0;
xr = R0 * sin(az); yr = R0 * cos(az);

Vrel_sq = du^2 + dv^2;
if Vrel_sq < 1e-6
    TCPA0 = inf; DCPA0 = R0;
else
    TCPA0 = -(xr*du + yr*dv) / Vrel_sq;
    DCPA0 = abs(xr*dv - yr*du) / sqrt(Vrel_sq);
end

%% 3. 原始危险判定输出
disp(' ');
disp('==================== 原始CPA计算结果 ====================');
fprintf('TCPA(最近会遇时间) = %.4f h ，折合 %.2f 分钟\n', TCPA0, TCPA0*60);
fprintf('DCPA(最近会遇距离) = %.4f 海里\n', DCPA0);
fprintf('设定安全距离阈值 = %.2f 海里\n', Dsafe);
disp('------------------------------------------------------');
risk_flag = 0;
if TCPA0 > 0 && DCPA0 < Dsafe
    risk_flag = 1;
    disp('⚠️ 危险判定：存在紧迫碰撞危险，启动会遇分析+遗传避让优化！');
elseif TCPA0 <= 0
    disp('✅ 危险判定：两船已驶过最近会遇点，持续远离，无需避让');
else
    disp('✅ 危险判定：会遇距离大于安全阈值，航行安全，无需避让');
end

%% 4. COLREGs会遇局面判定 & 确认让路船
give_way_flag = 0; % 1=本船让路；0=目标船让路；2=两船均让路
meet_type = '';
if risk_flag == 1
    disp(' ');
    disp('==================== 国际海上避碰规则 会遇局面分析 ====================');
    azi = azi0;
    if (azi >= 0 && azi <= 5) || (azi >= 355 && azi <= 360)
        meet_type = '对遇局面';
        give_way_flag = 2;
        disp('会遇局面：对遇局面');
        disp('让路责任：本船与目标船均为让路船，两船均需右转避让');
    elseif azi >= 112.5 && azi <= 247.5
        meet_type = '追越局面';
        give_way_flag = 1;
        disp('会遇局面：追越局面');
        disp('让路责任：本船为追越船，本船承担全部让路义务，右转避让');
    else
        if azi > 5 && azi < 112.5
            meet_type = '交叉相遇(目标在本船右舷)';
            give_way_flag = 1;
            disp('会遇局面：交叉相遇，目标船在本船右舷');
            disp('让路责任：本船为让路船，右转从目标船尾后方通过');
        elseif azi > 247.5 && azi < 355
            meet_type = '交叉相遇(目标在本船左舷)';
            give_way_flag = 0;
            disp('会遇局面：交叉相遇，目标船在本船左舷');
            disp('让路责任：目标船为让路船，本船保向保速，无需优化避让路线');
        end
    end
    disp('======================================================================');
end

%% 5. 仅当本船/两船需让路时，启动帕累托双目标遗传算法避让优化
BestX = [];
FitHist = zeros(100,2); % 预分配收敛记录数组，消除扩容警告
if risk_flag == 1 && give_way_flag ~= 0
    %% 5.1 决策变量上下限（遵循COLREG右转≥30°）
    Td_min = 0;    Td_max = TCPA0/2;    % 延迟直航时间不超过一半原始TCPA
    Ca_min = 30;   Ca_max = 70;         % 避让右转角度30~70°
    Ta_min = 0.05; Ta_max = TCPA0;      % 避让航行时间
    Cb_min = -70;  Cb_max = -30;        % 复航向左回正30~70°
    VarMin = [Td_min, Ca_min, Ta_min, Cb_min];
    VarMax = [Td_max, Ca_max, Ta_max, Cb_max];
    nVar = 4;        % 4个决策变量
    nPop = 50;       % 初代种群50
    MaxGen = 100;    % 迭代100次
    pCross = 0.7;    % 交叉概率
    pMut = 0.1;      % 变异概率
    EliteNum = 25;   % 轮盘筛选精英25
    OffNum = 25;     % 子代25，合并后种群维持50

    %% 5.2 初始化种群（预分配数组）
    Pop = zeros(nPop, nVar);
    for i = 1:nPop
        Pop(i,:) = VarMin + rand(1,nVar).*(VarMax - VarMin);
    end

    %% 5.3 GA主迭代循环
    for gen = 1:MaxGen
        % 计算种群双适应度（全部参数传参，无全局）
        Fit = zeros(nPop,2);
        for i = 1:nPop
            x = Pop(i,:);
            [f1,f2] = ObjFun(x, own_spd, own_crs, tar_spd, tar_crs, xr, yr, du, dv, deg2rad, TCPA0, Dsafe);
            Fit(i,:) = [f1,f2];
        end
        FitHist(gen,:) = mean(Fit,1);

        % 1) 轮盘赌选择精英25个
        ElitePop = RouletteSelect(Pop, Fit, EliteNum);
        % 2) 交叉变异生成25个子代（预分配子代数组）
        OffPop = zeros(OffNum, nVar);
        for k = 1:2:OffNum
            p1 = ElitePop(randi(EliteNum),:);
            p2 = ElitePop(randi(EliteNum),:);
            [c1,c2] = Cross(p1,p2,pCross,VarMin,VarMax);
            c1 = Mutate(c1,pMut,VarMin,VarMax);
            c2 = Mutate(c2,pMut,VarMin,VarMax);
            OffPop(k,:) = c1;
            if k+1 <= OffNum
                OffPop(k+1,:) = c2;
            end
        end
        % 3) 精英+子代合并新种群50
        Pop = [ElitePop; OffPop];
    end

    %% 5.4 最终种群帕累托非支配排序，提取最优折中解
    FinalFit = zeros(nPop,2);
    for i = 1:nPop
        FinalFit(i,:) = ObjFun(Pop(i,:), own_spd, own_crs, tar_spd, tar_crs, xr, yr, du, dv, deg2rad, TCPA0, Dsafe);
    end
    ParetoIdx = NonDominatedSort(FinalFit);
    ParetoPop = Pop(ParetoIdx,:);
    ParetoFit = FinalFit(ParetoIdx,:);
    % 折中最优解：距离原点欧氏距离最小的帕累托解
    Dist = sqrt(ParetoFit(:,1).^2 + ParetoFit(:,2).^2);
    [~,BestIdx] = min(Dist);
    BestX = ParetoPop(BestIdx,:);
    BestFit = ParetoFit(BestIdx,:);

    %% 5.5 输出最优避让方案参数
    disp(' ');
    disp('==================== 帕累托最优避让方案输出 ====================');
    fprintf('直航延迟时间Td = %.4f h (%.2f min)\n',BestX(1),BestX(1)*60);
    fprintf('避让右转角度Ca = %.2f °\n',BestX(2));
    fprintf('避让航行时间Ta = %.4f h (%.2f min)\n',BestX(3),BestX(3)*60);
    fprintf('复航回正角度Cb = %.2f °\n',BestX(4));
    fprintf('安全目标值(DCPA+TCPA) f1 = %.4f\n',BestFit(1));
    fprintf('经济目标值(总避让航程) f2 = %.4f n mile\n',BestFit(2));
    disp('================================================================');

    %% 6. 绘图1：迭代收敛曲线
    figure('Name','遗传算法迭代收敛曲线');
    subplot(2,1,1);
    plot(1:MaxGen,FitHist(:,1),'r-','LineWidth',1.2);
    xlabel('迭代次数');ylabel('平均安全目标f1=DCPA+TCPA');grid on;
    title('安全目标收敛曲线');
    subplot(2,1,2);
    plot(1:MaxGen,FitHist(:,2),'b-','LineWidth',1.2);
    xlabel('迭代次数');ylabel('平均经济目标f2=总避让航程');grid on;
    title('经济目标收敛曲线');

    %% 7. 绘图2：动态避让动画（避让船点 + 目标船安全警戒圆）
    DynamicAvoidAnim(BestX, own_spd, own_crs, tar_spd, tar_crs, xr, yr, deg2rad, TCPA0, Dsafe);
elseif risk_flag == 1 && give_way_flag == 0
    disp('无需启动避让优化：目标船为让路船，本船保向保速航行');
end

%% ====================== 子函数1：双目标适应度函数 ======================
function [f1,f2] = ObjFun(x, own_spd, own_crs, tar_spd, tar_crs, xr, yr, du, dv, deg2rad, TCPA0, Dsafe)
Td = x(1); Ca = x(2); Ta = x(3); Cb = x(4);
th0 = own_crs * deg2rad;

%% 阶段1：直航Td，本船原始航向
x0_1 = own_spd*sin(th0)*Td;
y0_1 = own_spd*cos(th0)*Td;
%% 阶段2：避让右转Ca，航向th0+Ca*deg2rad
th_a = th0 + Ca*deg2rad;
x0_2 = own_spd*sin(th_a)*Ta;
y0_2 = own_spd*cos(th_a)*Ta;
%% 阶段3：复航回正Cb，航向恢复原航向
th_b = th_a + Cb*deg2rad;
T_total = Td + Ta;
%% 延长仿真时长至完全驶过会遇点
T_sim = T_total + TCPA0 * 1.2;

%% 目标船全程匀速直线
xt_t = xr + tar_spd*sin(tar_crs*deg2rad)*T_sim;
yt_t = yr + tar_spd*cos(tar_crs*deg2rad)*T_sim;
%% 避让后本船终点坐标
x0_end = x0_1 + x0_2;
y0_end = y0_1 + y0_2;
%% 复航阶段继续航行至完全错开
T_back = TCPA0 * 1.2;
x0_back = own_spd*sin(th_b)*T_back;
y0_back = own_spd*cos(th_b)*T_back;
x0_final = x0_end + x0_back;
y0_final = y0_end + y0_back;

%% 计算避让后新DCPA/TCPA
dx = xt_t - x0_final; dy = yt_t - y0_final;
Vrel_sq = du^2 + dv^2;
TCPA_new = -(dx*du + dy*dv)/Vrel_sq;
DCPA_new = abs(dx*dv - dy*du)/sqrt(Vrel_sq);

%% 安全惩罚逻辑：不安全个体施加大惩罚，保证收敛下降
if TCPA_new < 0 || DCPA_new < Dsafe
    f1 = 200; % 不满足安全距离，高额惩罚
else
    f1 = DCPA_new + TCPA_new; % 安全目标：最小化DCPA+TCPA
end

%% 经济目标：直接最小化总绕行航程，曲线可正常收敛
L_total_avoid = own_spd * (Td + Ta + T_back);
f2 = L_total_avoid;
end

%% ====================== 子函数2：轮盘赌选择 ======================
function ElitePop = RouletteSelect(Pop,Fit,EliteNum)
% 防除零处理，避免计算报错
FitSum = Fit(:,1) + Fit(:,2) + 1e-6;
TotalFit = sum(1./FitSum); % 越小适应度越好，取倒数
Prob = (1./FitSum)/TotalFit;
ElitePop = zeros(EliteNum,size(Pop,2));
for i = 1:EliteNum
    r = rand();
    CumsumP = 0;
    for j = 1:size(Pop,1)
        CumsumP = CumsumP + Prob(j);
        if r <= CumsumP
            ElitePop(i,:) = Pop(j,:); break;
        end
    end
end
end

%% ====================== 子函数3：实数交叉 ======================
function [c1,c2] = Cross(p1,p2,pCross,VarMin,VarMax)
nVar = length(p1);
c1 = p1; c2 = p2;
if rand() < pCross
    alpha = rand();
    c1 = alpha*p1 + (1-alpha)*p2;
    c2 = (1-alpha)*p1 + alpha*p2;
end
% 边界约束
c1 = max(min(c1,VarMax),VarMin);
c2 = max(min(c2,VarMax),VarMin);
end

%% ====================== 子函数4：实数变异 ======================
function NewInd = Mutate(Ind,pMut,VarMin,VarMax)
nVar = length(Ind);
NewInd = Ind;
for i = 1:nVar
    if rand() < pMut
        NewInd(i) = VarMin(i) + rand()*(VarMax(i)-VarMin(i));
    end
end
end

%% ====================== 子函数5：快速非支配排序（帕累托） ======================
function ParetoIdx = NonDominatedSort(FitMat)
n = size(FitMat,1);
DomCount = zeros(n,1);
DomSet = cell(n,1);
ParetoIdx = [];
for p = 1:n
    for q = 1:n
        if p == q; continue; end
        % 判断p支配q：f1更小且f2更小
        if FitMat(p,1)<FitMat(q,1) && FitMat(p,2)<FitMat(q,2)
            DomSet{p} = [DomSet{p}, q];
        elseif FitMat(q,1)<FitMat(p,1) && FitMat(q,2)<FitMat(p,2)
            DomCount(p) = DomCount(p)+1;
        end
    end
    if DomCount(p) == 0
        ParetoIdx = [ParetoIdx, p];
    end
end
end

%% ====================== 子函数6：动态避让动画函数（核心新增） ======================
function DynamicAvoidAnim(BestX, own_spd, own_crs, tar_spd, tar_crs, xr, yr, deg2rad, TCPA0, Dsafe)
Td = BestX(1); Ca = BestX(2); Ta = BestX(3); Cb = BestX(4);
th0 = own_crs * deg2rad;
th_a = th0 + Ca*deg2rad;
th_b = th_a + Cb*deg2rad;
T_avoid = Td + Ta;
T_full_sim = T_avoid + TCPA0 * 1.5; % 完整仿真时长，覆盖错开全过程
t_step = linspace(0,T_full_sim,400); % 动画帧数

% 预分配坐标数组
x_own = zeros(size(t_step)); y_own = zeros(size(t_step));
x_tar = zeros(size(t_step)); y_tar = zeros(size(t_step));
tar_u = tar_spd*sin(tar_crs*deg2rad);
tar_v = tar_spd*cos(tar_crs*deg2rad);

% 目标船全程直线坐标
for k = 1:length(t_step)
    t = t_step(k);
    x_tar(k) = xr + tar_u * t;
    y_tar(k) = yr + tar_v * t;
end

% 避让船三段航行坐标
for k = 1:length(t_step)
    t = t_step(k);
    if t <= Td
        x_own(k) = own_spd*sin(th0)*t;
        y_own(k) = own_spd*cos(th0)*t;
    elseif t > Td && t <= T_avoid
        dt = t - Td;
        x_own(k) = own_spd*sin(th0)*Td + own_spd*sin(th_a)*dt;
        y_own(k) = own_spd*cos(th0)*Td + own_spd*cos(th_a)*dt;
    else
        dt1 = Ta;
        dt2 = t - T_avoid;
        x_own(k) = own_spd*sin(th0)*Td + own_spd*sin(th_a)*dt1 + own_spd*sin(th_b)*dt2;
        y_own(k) = own_spd*cos(th0)*Td + own_spd*cos(th_a)*dt1 + own_spd*cos(th_b)*dt2;
    end
end

%% 新建动画窗口
figure('Name','船舶避让动态动画');
hold on; grid on; axis equal;
xlabel('东向距离(n mile)');ylabel('北向距离(n mile)');
title('动态避让动画：蓝色点=避让船 | 红色点+虚线圆=目标船（安全警戒圈）');

% 绘制完整历史航迹（静态背景）
plot(x_own,y_own,'b-','LineWidth',1,'DisplayName','避让船完整航迹');
plot(x_tar,y_tar,'r-','LineWidth',1,'DisplayName','目标船直航迹');

% 动态对象初始化
h_own_point = plot(0,0,'bo','MarkerSize',9,'MarkerFaceColor','b','DisplayName','避让船');
h_tar_point = plot(xr,yr,'ro','MarkerSize',9,'MarkerFaceColor','r','DisplayName','目标船');
% 安全警戒圆（360度采样）
theta_circle = linspace(0,2*pi,100);
circ_x = Dsafe * cos(theta_circle);
circ_y = Dsafe * sin(theta_circle);
h_safe_circle = plot(xr+circ_x, yr+circ_y, 'r--','LineWidth',1.2,'DisplayName','安全警戒圈');
legend('Location','best');

%% 逐帧动画循环
for frame = 1:length(t_step)
    % 更新避让船点位
    set(h_own_point, 'XData',x_own(frame), 'YData',y_own(frame));
    % 更新目标船点位与警戒圆
    tx = x_tar(frame); ty = y_tar(frame);
    set(h_tar_point, 'XData',tx, 'YData',ty);
    set(h_safe_circle, 'XData',tx + circ_x, 'YData',ty + circ_y);
    pause(0.02); % 控制动画速度，数值越小播放越快
end
hold off;
end