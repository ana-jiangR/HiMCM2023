clc;clear;close all;
%% 数据提取
data_ori=readtable('Busline_Stops.xlsx',"TextType","string","TreatAsMissing",[".","NA"]);
ind_missing=ismissing(data_ori);%缺省数据序号
B=[];
for i=1:size(data_ori,1)
    line_ori(i).lineid=data_ori{i,1};
    line_ori(i).distance=data_ori{i,2};
    A=[false false false ~ind_missing(i,4:end)];
    line_ori(i).station=data_ori{i,A};
    B=[B;line_ori(i).station'];%所有站点
end
C=unique(B);%站点名称及编号
%% 依据站点编号使线路数字化
for i=1:size(line_ori,2)
    [is,pos]=ismember(line_ori(i).station,C);
    line_ori(i).ind_station=pos;
end

%% 参数设置
Capacity=250;%容量，Kwh
Fuel_Consumption=0.65;%燃料消耗速度，kwh/km
Charge_speed=60;%充电速度kwh/h
T_layover=5;%单次停靠时间，为5minute
C_vehicle= 4196602;%车辆价格sek
C_Electric=0.172;%电量价格sek/km
C_Diesel=5.714;%燃油价格sek/km
C_infrastructure=2115000;%充电站sek/statuon
SOC=0.3;%最低能量状态
V_vehicle=40;%公交速度25km/h
H=1/3;%班次间隔,h
fund=4e8;%公交系统营业收入
Charge_kwh=Charge_speed*T_layover/60;%单次停靠充电

%% 排除掉绝对不可能的路线，计算每条路线的需要投放的电动公交车数量
ind_impossible=[];%保存不可能的序号
Station=zeros(size(line_ori,2),length(C));
for i=1:size(line_ori,2)
    i
    A=[];
    B=[];
    ind_station=line_ori(i).ind_station;
    col=length(ind_station);
    d=line_ori(i).distance/col;
    A=zeros(2*col,col);
    for j=1:2*col
        if j<=col
            A(j,1:j)=1;
        else
            A(j,:)=A(j-1,:);
            A(j,2*col-j+1)=2;
        end
    end
    %不等式约束系数，AX+C-D《C；AX+C-D》soc*C;
    A=[A;-A]*Charge_kwh;%系数
    B=(1:2*col)'*d*Fuel_Consumption;
    B=[B;Capacity-SOC*Capacity-(1:2*col)'*d*Fuel_Consumption];
    f=ones(1,col);
    intcon=1:col;
    Aeq=[];
    Beq=[];
    lb=zeros(1,col);
    ub=ones(1,col);
    [x,fval]=intlinprog(f,intcon,A,B,Aeq,Beq,lb,ub);
    if isempty(x)
        ind_impossible=[ind_impossible i];
        X(i)=0;%不可行的设置为0
        continue
    end
    %找到实际对应的站点
    Station(i,ind_station(logical(x)))=1;

    %每条线路多少车，至少满足每20min一班车H-2*D/V+H*(x-1)>(2D*cons-Ax)/Vcharge
    intcon=1;
    B=-2*line_ori(i).distance/V_vehicle+T_layover/60*sum(x)+ ...
        -(2*line_ori(i).distance*Fuel_Consumption-A(end,:)*x*Charge_kwh)/Charge_speed;
    x=[];
    A=-H;
    Aeq=[];
    Beq=[];
    lb=[];
    ub=[];
    [X(i),fval]=intlinprog(1,intcon,A,B,Aeq,Beq,lb,ub);
end
line_num_vehicle=X;%每条路线的电动公交车数
total_num_vehicle=sum(X);%电动公交车总数
line_sum_Estaton=sum(Station,2)'+line_num_vehicle;%认为线路多少辆车则在基站需配备多少充电设置

Final_line_station_Total=[];
Final_line_vehicle=[];
Final_charger_station_Total=[];
Final_Line=[];

%% 考虑成本对线路进行规划，目标函数为最大化电动公交线路里程（里程越长，电动车越优）,线路为01矩阵
f=-data_ori{:,2};%线路里程,标准化，取负
intcon=1:length(f);
%车辆价格+充电站价格（有重复）+燃料节省价格
A=line_num_vehicle*C_vehicle+line_sum_Estaton*C_infrastructure+365*10/H*(C_Electric-C_Diesel)*data_ori{:,2}';
A0=line_num_vehicle*C_vehicle+line_sum_Estaton*C_infrastructure+365*10/H;
ind_A0=ones(1,length(line_sum_Estaton));
ind_A0(ind_impossible)=0;
C_all=A0*ind_A0'*0.5;%全部转化的资金
A(ind_impossible)=1e10;%使不可能线路不可规划
B=C_all+fund;
Aeq=[];
Beq=[];
lb=zeros(1,length(f));
ub=ones(1,length(f));
[line_allow,fval]=intlinprog(f,intcon,A,B,Aeq,Beq,lb,ub);%求解被规划的线路
fval=-fval;%线路总长

%% 确定最终的站点规划和总资金
%确定数据
final_Line=line_allow;%最终所确定的线路
final_line_station_Total=line_sum_Estaton;%各线路所需充电设施数，基地+路线中
final_line_station_Total(~logical(line_allow))=0;%未规划的线路
final_line_vehicle=line_num_vehicle;%各线路所需电动汽车数
final_line_vehicle(~logical(line_allow))=0;%未规划的线路
final_charger_station_Total=sum(Station(logical(line_allow),:));%各个站点需要安装的充电设施数量
%确定资金
final_cost=A*line_allow;%应为2*fund以下
Final_line_station_Total=[Final_line_station_Total;final_line_station_Total];
Final_line_vehicle=[Final_line_vehicle;final_line_vehicle];
Final_charger_station_Total=[Final_charger_station_Total;final_charger_station_Total];
Final_Line=[Final_Line final_Line];

figure






















