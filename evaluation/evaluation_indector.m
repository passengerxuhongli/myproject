
clear
clc
%%
% Trj_real = xlsread('link-ebike.xlsx',1);
% for i= 1:length(index)-1
%    Trj_real(find(Trj_real(:,1)==Trj_real(index(i),1)):find(Trj_real(:,1)==Trj_real(index(i+1),1))-1,1) = Trj_real(index(i),1);
% end%把轨迹编号补全
% Trj_real(isnan(Trj_real(:,2)),:)=[];%剔除标题行
% t0 = min(Trj_real(:,2));%仿真初始时刻
% Trj_real(:,2) = (Trj_real(:,2)-t0)*24*3600;%同一时间刻度为0.12秒的倍数
% % 把link-ebike.xlsx数据格式变成‘评价基础数据’的格式
%读取数据
Trj_real = xlsread('评价基础数据.xlsx',1);
Trj_sfm = xlsread('评价基础数据.xlsx',2);
Trj_MSFMcomfortzone = xlsread('评价基础数据.xlsx',3);
index = find(~isnan(Trj_real(:,1)));%real
%% ================ 超车效果检验 ==============================
%% 超车次数，超车横纵向间距等参数提取 
%field data 超车信息
overtake_number=0;bike_number=[];left_field = 0;times = 0;%左侧超车 
time_overtake = []; longitudinal_distance =[]; lateral_space =[];lateral_distance =[];tip =0 ;%超车记录标签
for i=1:max(Trj_real(:,1))%后车车辆循环
    for T= min(Trj_real(Trj_real(:,1)==i,2)):0.12:max(Trj_real(Trj_real(:,1)==i,2))-0.12%单辆车运行时间遍历
        front_bike = Trj_real(abs(Trj_real(:,2)-T)<0.06,:);%寻找所有前车         
        back_bike = Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-T)<0.06,:);%寻找当前后车
        front_bike(front_bike(:,1)==i,: )=[]; %剔除本身车辆        
        front_bike1 = Trj_real(abs(Trj_real(:,2)-T-0.12)<0.06,:);%寻找下一步长所有前车
        back_bike1 = Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-T-0.12)<0.06,:);%寻找下一步长的后车
        front_bike1(front_bike1(:,1)==i,: )=[]; %剔除本身车辆
        if ~isempty(front_bike(:,4))
            for j=1:length(front_bike(:,1))
             if  ~isempty( find(front_bike1(:,1)==front_bike(j,1)))
                if  (back_bike(1,4)-front_bike(j,4))<0 &(back_bike1(1,4)-front_bike1(front_bike1(:,1)==front_bike(j,1),4))>=0  %超车？  %判断该车辆是否为前车&是否发生后车超越前车现象？
                    for t=max(min(Trj_real(Trj_real(:,1)==i,2)),min(Trj_real(Trj_real(:,1)==front_bike(j,1),2))):0.12:T
                        zone = Trj_real(Trj_real(:,4)>=Trj_real(Trj_real(:,1)==front_bike(j,1)&abs(Trj_real(:,2)-t)<0.06,4)-2.5&Trj_real(:,4)<=Trj_real(Trj_real(:,1)==front_bike(j,1)&abs(Trj_real(:,2)-t)<0.06,4)+2.5,:);%统计区域为5m
                        density = length(zone(abs(zone(:,2)-t)<0.06,1))/(4*5);
                        r1 = abs( 15.3911+0.5*Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,6)+0.1784*Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,6).^2-0.4234*Trj_real(Trj_real(:,1)==front_bike(j,1)&abs(Trj_real(:,2)-t)<0.06,6).^2-56.1606*density ); %长轴
                        zone1 = Trj_real(Trj_real(:,4)>=Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,4)-5&Trj_real(:,4)<=Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,4)+5,:);%统计区域为10m
                        density1 = length(zone1(abs(zone1(:,2)-t)<0.06,1))/(4*10);
                        interdistance = min(0.156*density1^(-0.18746),1.6);
                        comfort_lateral = min(0.8+interdistance,1.6);
                        r2 = 0.8-0.5*Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,6)*sin(pi-2*acos(9.*1.0*tan(pi*15/180)/(Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,6).^2)))+comfort_lateral;
                        a1=Trj_real(abs(Trj_real(:,2)-t)<0.06&Trj_real(:,1)==i,4);
                        a2=Trj_real(abs(Trj_real(:,2)-t)<0.06&Trj_real(:,1)==front_bike(j,1),4);
                        if abs(a1(1)-a2(1))<=r1
                            a3=Trj_real(abs(Trj_real(:,2)-t)<0.06&Trj_real(:,1)==i,3);
                            a4=Trj_real(abs(Trj_real(:,2)-t)<0.06&Trj_real(:,1)==front_bike(j,1),3);
                            if abs(a3(1)-a4(1))<=r2%侵犯舒适空间？                               
%                                delta_x = abs(front_bike(j,3)-back_bike(1,3));
%                                delta_y = abs(front_bike(j,4)-back_bike(1,4));
%                                 yy = delta_y.^2/r1 + delta_x.^2/r2; 
%                                 if yy<=1 & (front_bike(j,4)-back_bike(1,4)>=0)
                                overtake_number = overtake_number+1;%超车次数累加
                                time_overtake =[time_overtake; t];%超车时刻
                                longitudinal_distance =[ longitudinal_distance; abs(a2(1)-a1(1))];%纵向间距m
                                lateral_space =[ lateral_space; max(a4(1)-2,6-a4(1))];%超车开始时前车与边界之间的间距最大值
                                lateral_distance =[lateral_distance; abs( front_bike1(j,3)-back_bike1(1,3))];%超车并行点处，两车间距。
                                if front_bike(j,3)-back_bike(1,3)>=0
                                    left_field = left_field+1;
                                end
                                tip = 1;
                                bike_number=[bike_number;i];%超车车辆-后车
                                break;
%                                 end 
                            end
                        end
                    end
                    times =times+1;
                end
            end
                if tip==1
                    break;
                    tip=0;
                end
            end
        end
    end
end
    %输出参数
    times
    overtake_number
    left_field
    xlswrite('超车参数1',[bike_number,time_overtake,longitudinal_distance,lateral_space,lateral_distance],1);

    %sfm 超车次数等超车信息
    overtake_number=0;bike_number=[];left_sfm = 0;times = 0;%左侧超车
    time_overtake = []; longitudinal_distance =[]; lateral_space =[];lateral_distance =[];tip =0 ;%超车记录标签
    for i=1:max(Trj_sfm(:,1))%后车车辆循环
        for T= min(Trj_sfm(Trj_sfm(:,1)==i,2)):0.12:max(Trj_sfm(Trj_sfm(:,1)==i,2))-0.12%单辆车运行时间遍历
            front_bike = Trj_sfm(abs(Trj_sfm(:,2)-T)<0.06,:);%寻找所有前车
            back_bike = Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-T)<0.06,:);%寻找当前后车
            front_bike(front_bike(:,1)==i,: )=[]; %剔除本身车辆
            front_bike1 = Trj_sfm(abs(Trj_sfm(:,2)-T-0.12)<0.06,:);%寻找下一步长所有前车
            back_bike1 = Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-T-0.12)<0.06,:);%寻找下一步长的后车
            front_bike1(front_bike1(:,1)==i,: )=[]; %剔除本身车辆
            if ~isempty(front_bike(:,4))
                for j=1:length(front_bike(:,1))
                 if  ~isempty( find(front_bike1(:,1)==front_bike(j,1)))
                    if   (back_bike(1,4)-front_bike(j,4))<0 & (back_bike1(1,4)-front_bike1(front_bike1(:,1)==front_bike(j,1),4))>=0  %超车？
                        for t= max(min(Trj_sfm(Trj_sfm(:,1)==i,2)),min(Trj_sfm(Trj_sfm(:,1)==front_bike(j,1),2))):0.12:T                                                       
                        zone = Trj_sfm(Trj_sfm(:,4)>=Trj_sfm(Trj_sfm(:,1)==front_bike(j,1)&abs(Trj_sfm(:,2)-t)<0.06,4)-2.5&Trj_sfm(:,4)<=Trj_sfm(Trj_sfm(:,1)==front_bike(j,1)&abs(Trj_sfm(:,2)-t)<0.06,4)+2.5,:);%统计区域为5m
                        density = length(zone(abs(zone(:,2)-t)<0.06,1))/(4*5);
                        r1 = abs( 15.3911+0.5*Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,6)+0.1784*Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,6).^2-0.4234*Trj_sfm(Trj_sfm(:,1)==front_bike(j,1)&abs(Trj_sfm(:,2)-t)<0.06,6).^2-56.1606*density ); %长轴
                        zone1 = Trj_sfm(Trj_sfm(:,4)>=Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,4)-5&Trj_sfm(:,4)<=Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,4)+5,:);%统计区域为10m
                        density1 = length(zone1(abs(zone1(:,2)-t)<0.06,1))/(4*10);
                        interdistance = min(0.156*density1^(-0.18746),1.6);
                        comfort_lateral = min(0.6+interdistance,1.6);
                        r2 = 0.8-0.5*Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,6)*sin(pi-2*acos(9.*1.0*tan(pi*15/180)/(Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,6).^2)))+comfort_lateral;                            
                            a1=Trj_sfm(abs(Trj_sfm(:,2)-t)<0.06&Trj_sfm(:,1)==i,4);
                            a2=Trj_sfm(abs(Trj_sfm(:,2)-t)<0.06&Trj_sfm(:,1)==front_bike(j,1),4);
                            if abs(a1(1)-a2(1))<=r1
                                a3=Trj_sfm(abs(Trj_sfm(:,2)-t)<0.06&Trj_sfm(:,1)==i,3);
                                a4=Trj_sfm(abs(Trj_sfm(:,2)-t)<0.06&Trj_sfm(:,1)==front_bike(j,1),3);
                                if abs(a3(1)-a4(1))<=r2%侵犯舒适空间？
                                    
                                    delta_x = abs(front_bike(j,3)-back_bike(1,3));
                                    delta_y = abs(front_bike(j,4)-back_bike(1,4));
                                    yy = delta_y.^2/r1 + delta_x.^2/r2; 
                                    if yy<=1 & (front_bike(j,4)-back_bike(1,4)>=0)                                    
                                    overtake_number = overtake_number+1;%超车次数累加
                                    time_overtake =[time_overtake; t];%超车时刻
                                    longitudinal_distance =[ longitudinal_distance; abs(a2(1)-a1(1))];%纵向间距m
                                    lateral_space =[ lateral_space; max(a4(1)-2,6-a4(1))];%超车开始时前车与边界之间的间距最大值
                                    lateral_distance =[lateral_distance; abs( front_bike1(j,3)-back_bike1(1,3))];%超车并行点处，两车间距。
                                    if front_bike(j,3)-back_bike(1,3)>=0
                                        left_sfm = left_sfm+1;
                                    end
                                    tip = 1;
                                    bike_number=[bike_number;i];%超车车辆-后车
                                    break;
                                    end 
                                end
                            end
                        end
                        times = times+1;
                    end
                end
                    if tip==1
                        break;
                        tip=0;
                    end
                end
            end
        end
    end
    %输出参数
    times
    overtake_number
    left_sfm
    xlswrite('超车参数1',[bike_number,time_overtake,longitudinal_distance,lateral_space,lateral_distance],2);

    %msfmcomfortzone 超车次数等超车信息
    overtake_number=0;bike_number=[];left_msfm = 0;times = 0;%左侧超车
    time_overtake = []; longitudinal_distance =[]; lateral_space =[];lateral_distance =[];tip =0 ;%超车记录标签
    for i=1:max(Trj_MSFMcomfortzone(:,1))%后车车辆循环
        for T= min(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,2)):0.12:max(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,2))-0.12%单辆车运行时间遍历
            
            front_bike = Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-T)<0.06,:);%寻找所有前车
            back_bike = Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-T)<0.06,:);%寻找当前后车
            front_bike(front_bike(:,1)==i,: )=[]; %剔除本身车辆
            front_bike1 = Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-T-0.12)<0.06,:);%寻找下一步长所有前车
            back_bike1 = Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-T-0.12)<0.06,:);%寻找下一步长的后车
            front_bike1(front_bike1(:,1)==i,: )=[]; %剔除本身车辆

            if ~isempty(front_bike(:,4))
                for j=1:length(front_bike(:,1))
                  if  ~isempty( find(front_bike1(:,1)==front_bike(j,1)))
                    if   (back_bike(1,4)-front_bike(j,4))<0 & (back_bike1(1,4)-front_bike1(front_bike1(:,1)==front_bike(j,1),4))>=0  %超车？
                        for t= max(min(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,2)),min(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==front_bike(j,1),2))):0.12:T   
                        zone = Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,4)>=Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==front_bike(j,1)...
                            &abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,4)-2.5&Trj_MSFMcomfortzone(:,4)<=Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==front_bike(j,1)...
                            &abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,4)+2.5,:);%统计区域为5m
                        density = length(zone(abs(zone(:,2)-t)<0.06,1))/(4*5);
                        r1 = abs( 15.3911+0.5*Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6)+0.1784*Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6).^2-0.4234*Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==front_bike(j,1)&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6).^2-56.1606*density ); %长轴
                        zone1 = Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,4)>=Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,4)-5&Trj_MSFMcomfortzone(:,4)<=Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,4)+5,:);%统计区域为10m
                        density1 = length(zone1(abs(zone1(:,2)-t)<0.06,1))/(4*10);
                        interdistance = min(0.156*density1^(-0.18746),1.6);
                        comfort_lateral = min(0.6+interdistance,1.6);
                        r2 = 0.8-0.5*Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6)*sin(pi-2*acos(9.*1.0*tan(pi*15/180)/(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6).^2)))+comfort_lateral;                            
                            a1=Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-t)<0.06&Trj_MSFMcomfortzone(:,1)==i,4);
                            a2=Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-t)<0.06&Trj_MSFMcomfortzone(:,1)==front_bike(j,1),4);
                            if abs(a1(1)-a2(1))<=r1
                                a3=Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-t)<0.06&Trj_MSFMcomfortzone(:,1)==i,3);
                                a4=Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-t)<0.06&Trj_MSFMcomfortzone(:,1)==front_bike(j,1),3);
                                if abs(a3(1)-a4(1))<=r2%侵犯舒适空间？
                                    delta_x = abs(front_bike(j,3)-back_bike(1,3));
                                    delta_y = abs(front_bike(j,4)-back_bike(1,4));
                                    yy = delta_y.^2/r1 + delta_x.^2/r2; 
                                    if yy<=1 & (front_bike(j,4)-back_bike(1,4)>=0)                                    
                                    overtake_number = overtake_number+1;%超车次数累加
                                    time_overtake =[time_overtake; t];%超车时刻
                                    longitudinal_distance =[ longitudinal_distance; abs(a2(1)-a1(1))];%纵向间距m
                                    lateral_space =[ lateral_space; max(a4(1)-2,6-a4(1))];%超车开始时前车与边界之间的间距最大值
                                    lateral_distance =[lateral_distance; abs( front_bike1(j,3)-back_bike1(1,3))];%
                                    if  front_bike(j,3)-back_bike(1,3)>=0
                                        left_msfm = left_msfm+1;
                                    end
                                    tip = 1;
                                    bike_number=[bike_number;i];%超车车辆-后车
                                    break;
                                    end 
                                end
                            end
                        end
                        times = times+1;
                    end
                 end
                    if tip==1
                        break;
                        tip=0;
                    end
                end
            end
        end
    end
    %输出参数
    times
    overtake_number
    left_msfm
    xlswrite('超车参数1',[bike_number,time_overtake,longitudinal_distance,lateral_space,lateral_distance],3);

    %% =======================超车横纵向间距分布对比==================
    field_overtake = xlsread('超车参数1.xls',1);
    sfm_overtake = xlsread('超车参数1.xls',2);
    msfm_overtake = xlsread('超车参数1.xls',3);
    % 超车信息矩阵箱型图
    x1 = field_overtake(:,2:5);
    x2 = sfm_overtake(:,2:5);
    x3 = msfm_overtake(:,2:5);
    %x1补齐维度，与x2一致
    xx1 = zeros(length(x2),4);
    xx1(:,:) = NaN;
    xx1(1:length(x1),:) = x1;
    %x3补齐维度，与x2一致
    xx3 = zeros(length(x2),4);
    xx3(:,:) = NaN;
    xx3(1:length(x3),:) = x3;
    X = [xx1 x2 xx3];%综合三种情况下超车信息矩阵
%  % 超车参数箱型图
%  LY_FIELD = X(:,2);
%  LY_MYMODEL = X(:,10);
%  CLEARENCE_FIELD = X(:,3);
%  CLEATENCE_MYMODEL = X(:,11);   
%  figure
%     subplot(1,2,1)
%     boxplot([LY_FIELD  LY_MYMODEL],{'LY_FIELD', 'LY_MYMODEL'})
%     grid on;
%     subplot(1,2,2)
%     boxplot([CLEARENCE_FIELD CLEATENCE_MYMODEL],{'CLEAR_FIELD', 'CLEAR_MYMODEL'})
%     grid on;    
% %        saveas(gcf,['超车过程时刻箱型图','.fig']);    
%  LX_FIELD = X(:,4);
%  LX_MYMODEL = X(:,12)+0.8; 
% figure
%     boxplot([LX_FIELD LX_MYMODEL],{'LX_FIELD', 'LX_MYMODEL'})
%     grid on; 
% %        saveas(gcf,['超车并行时刻箱型图','.fig']);

%% 超车时刻箱型图
    figure
    subplot(1,2,1)
    boxplot([X(:,2) X(:,6) X(:,10)],{'FIELDlongitidual', 'SFMlongitidual','MYMODELlongitidual'})
    grid on;
    subplot(1,2,2)
    boxplot([X(:,3) X(:,7) X(:,11)],{'FIELDclearance', 'SFMclearance','MYMODELclearance'})
    grid on;
    saveas(gcf,['超车过程箱型图','.fig']);   
    
figure
    boxplot([X(:,4) X(:,8) X(:,12)+0.8],{'FIELDlateral', 'SFMlateral', 'MYMODELlateral'})
    grid on;
   saveas(gcf,['超车并行时刻箱型图','.fig']);

    
%% --------------------指标3瞬时位置均方误差，数据------------------------
% 整体横坐标标准误差：各测量值误差平方和的平均值求平方根。
    %MYMODEL
    for i = 1:max(Trj_real(:,1))
        cd = min([ sum(Trj_MSFMcomfortzone(:,1)==i),sum(Trj_real(:,1)==i)]);
        data_sfm_fix_passing = find(Trj_MSFMcomfortzone(:,1)==i);
        if cd<=0
            continue;
        end
        px_sfm_fix_passing(data_sfm_fix_passing(1):data_sfm_fix_passing(1)+cd-1,1) = (Trj_real(index(i)+1:index(i)+cd,3) - Trj_MSFMcomfortzone( data_sfm_fix_passing(1):data_sfm_fix_passing(1)+cd-1,3)).^2;
    end
    px_sfm_fix_passing = px_sfm_fix_passing(~isnan(px_sfm_fix_passing));
    error_Px_sfm_fix_passing = sqrt(sum(px_sfm_fix_passing)/length(px_sfm_fix_passing));%sfm_fix_passing
    %MYMODEL
    for i = 1:max(Trj_real(:,1))
        cd = min([ sum(Trj_MSFMcomfortzone(:,1)==i),sum(Trj_real(:,1)==i)]);
        data_sfm_fix_passing = find(Trj_MSFMcomfortzone(:,1)==i);
        if cd<=0
            continue;
        end
        py_sfm_fix_passing(data_sfm_fix_passing(1):data_sfm_fix_passing(1)+cd-1,1) = (Trj_real(index(i)+1:index(i)+cd,4) - Trj_MSFMcomfortzone( data_sfm_fix_passing(1):data_sfm_fix_passing(1)+cd-1,4)).^2;
    end
    py_sfm_fix_passing = py_sfm_fix_passing(~isnan(py_sfm_fix_passing));
    error_py_sfm_fix_passing = sqrt(sum(py_sfm_fix_passing)/length(py_sfm_fix_passing));%sfm_fix_passing
    
      %% 整体坐标均方误差  
        %sfm
    for i = 1:max(Trj_real(:,1))
        cd = min([ sum(Trj_sfm(:,1)==i),sum(Trj_real(:,1)==i)]);
        data_sfm = find(Trj_sfm(:,1)==i);
        if cd<=0
            continue;
        end
        px_sfm(data_sfm(1):data_sfm(1)+cd-1,1) = (Trj_real(index(i)+1:index(i)+cd,3) - Trj_sfm( data_sfm(1):data_sfm(1)+cd-1,3)).^2;
    end
    px_sfm = px_sfm(~isnan(px_sfm));
    error_Px_sfm = sqrt(sum(px_sfm)/length(px_sfm));%sfm
    %sfm
    for i = 1:max(Trj_real(:,1))
        cd = min([ sum(Trj_sfm(:,1)==i),sum(Trj_real(:,1)==i)]);
        data_sfm = find(Trj_sfm(:,1)==i);
        if cd<=0
            continue;
        end
        py_sfm(data_sfm(1):data_sfm(1)+cd-1,1) = (Trj_real(index(i)+1:index(i)+cd,4) - Trj_sfm( data_sfm(1):data_sfm(1)+cd-1,4)).^2;
    end
    py_sfm = py_sfm(~isnan(py_sfm));
    error_py = sqrt(sum(py_sfm)/length(py_sfm));%sfm_fix_passing    
    %% =====================超车效果图============================
    data = xlsread('link-ebike','受阻超车车辆编号测试集');
    figure;
    for i = 1:size(data,1)
        plot( Trj_real( Trj_real(:,1)==data(i,1),4),Trj_real( Trj_real(:,1)==data(i,1),3),'-b');axis([-10,90,1.5,7]);
        hold on;
        plot( Trj_sfm( Trj_sfm(:,1)==data(i,1),4),Trj_sfm( Trj_sfm(:,1)==data(i,1),3),'-g');axis([-10,90,1.5,7])
        plot(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==data(i,1),4),Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==data(i,1),3),'-r');axis([-10,90,1.5,7]);
        grid;
         h = legend('FIELD','SFM','MSFM');
%            h = legend('FIELD','MYMODEL');
        set(h, 'Box', 'off');
        set(h,'Interpreter','none')
        axis([-10 90 0 6]);
        title('overtaking-Sx-Sy');ylabel('lateral position /m');xlabel('longitudinal position /m');
%        saveas(gcf,['超车单个轨迹',num2str(i),'.fig']);
  saveas(gcf,['超车单个轨迹',num2str(i),'.png']);
             hold off;
    end    
    
    %% ======================轨迹位置图Sx-Sy ======================
    % 编号、时间、横向位置、纵向位置、横向仿真速度、纵向仿真速度、横向仿真加速度、纵向仿真加速度
    %--------------------指标1轨迹图------------------------
    %% 动态单条
    figure
    for i=193:401
        plot( Trj_real( Trj_real(:,1)==i,4),Trj_real( Trj_real(:,1)==i,3),'-b.');axis([-10,90,0,6]);
        hold on;
        plot( Trj_sfm( Trj_sfm(:,1)==i,4),Trj_sfm( Trj_sfm(:,1)==i,3),'-g.');axis([-10,90,1.5,7])
        plot(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,4),Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,3),'-r.');axis([-10,90,1.5,7]);
        grid;
        h = legend('FIELD','SFM','MYMODEL');
%         h = legend('FIELD','MYMODEL');
        set(h, 'Box', 'off');
        set(h,'Interpreter','none')
        axis([-10 90 0 6]);
        title('Sx-Sy');ylabel('lateral distance /m');xlabel('longitudinal distance /m');
        pause(0.5)
  hold off;
%         saveas(gcf,['单个轨迹',num2str(i),'.fig']);
 saveas(gcf,['单个轨迹',num2str(i),'.png']);
    end
    %% 总体轨迹
    figure
    hold on;
    for i =1:max(Trj_real(:,1))
        plot( Trj_real( Trj_real(:,1)==i,4),Trj_real( Trj_real(:,1)==i,3),'b');
        hold on;
        plot(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,4),Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,3),'r');
      
    end
    h = legend('FIELD','MYMODEL');
    set(h, 'Box', 'off');
    set(h,'Interpreter','none')
    axis([-10 90 1 8]);
    grid;
    title('Sx-Sy');ylabel('lateral distance /m');xlabel('longitudinal distance /m');
     saveas(gcf, '整体轨迹','fig');
    %% 速度指标
    %% --------------------区间速度分布-----------------------
    %% 区间平均速度分布，距离除以时间
    speed = zeros(length(Trj_real(:,1)),4);
    [C, ia, ic] = unique(Trj_real(:,1),'rows');%取ia
    for i = 1:size(ia,1)
        Trj_r = Trj_real(Trj_real(:,1)==Trj_real(ia(i),1),:);
        Trj_s = Trj_sfm(Trj_sfm(:,1)==Trj_real(ia(i),1),:);
        Trj_s_f_p = Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==Trj_real(ia(i),1),:);
        v_real = (Trj_r(end,4)-Trj_r(1,4)) / ((Trj_r(end,2)-Trj_r(1,2)));%real
        v_sfm =  (Trj_s(end,4)-Trj_s(1,4)) / (Trj_s(end,2)-Trj_s(1,2));%SFM
        v_sfm_fix_passing = (Trj_s_f_p(end,4)-Trj_s_f_p(1,4)) / (Trj_s_f_p(end,2)-Trj_s_f_p(1,2));%SFM_fix+passing
        speed(i,:) = [Trj_real(ia(i),1),v_real,v_sfm,v_sfm_fix_passing];
    end
    speed = speed(speed(:,1)~=0,:);
    
    V_real = fitdist(speed(:,2),'norm');
    V_sfm = fitdist(speed(:,3),'norm');
    V_sfm_fix_passing = fitdist(speed(:,4),'norm');
    
    figure;
    hold on
    x = [0:0.01:14];
    y1 = pdf(V_real,x);
    y2 = pdf(V_sfm,x);
    y3 = pdf(V_sfm_fix_passing,x);
    plot(x,y1);
    plot(x,y2);
    plot(x,y3);
    hold on
    [a1,b1]=hist(speed(:,2));%速度,,,作图后手动在工具箱中调整外观美化
    bar(b1,a1/sum(a1));
    [a2,b2]=hist(speed(:,3));%速度
    bar(b2,a2/sum(a2));
    [a3,b3]=hist(speed(:,4));%速度
    bar(b3,a3/sum(a3));
    title('The distribution of travel speed');ylabel('Frequency');xlabel('travel speed m/s');
       saveas(gcf, 'The distribution of travel speed','fig');
    %% ----------
    %% 每仿真秒距离 compute
    num = min([length(Trj_real),length(Trj_sfm),length(Trj_MSFMcomfortzone)]);
    t_start = Trj_MSFMcomfortzone(1,2);
    t_end = max(Trj_MSFMcomfortzone(:,2));
    L_REAL = zeros(num,2);%间距、平均速度
    L_SFM = zeros(num,2);
    L_SFM_FIX_PASSING = zeros(num,2);
    c_real =1;
    c_sfm = 1;
    c_sfm_fix_passing = 1;
    for t = t_start:0.12:t_end
        p_index_real = find(abs((Trj_real(:,2))-t)<0.06);
        p_index_sfm = find(abs(Trj_sfm(:,2)-t)<0.06);
        p_index_sfm_fix_passing = find(abs(Trj_MSFMcomfortzone(:,2)-t)<0.06);
        if length(p_index_sfm)>1
            if length(p_index_real)>1
                if length(p_index_sfm_fix_passing)>1
                    inter_length_real = sqrt(diff(Trj_real(p_index_real,3)).^2+diff(Trj_real(p_index_real,4)).^2);%两车间距%real
                    inter_length_sfm = sqrt(diff(Trj_sfm(p_index_sfm,3)).^2+diff(Trj_sfm(p_index_sfm,4)).^2);%两车间距%SFM
                    inter_length_sfm_fix_passing = sqrt(diff(Trj_MSFMcomfortzone(p_index_sfm_fix_passing,3)).^2+diff(Trj_MSFMcomfortzone(p_index_sfm_fix_passing,4)).^2);%两车间距%修改SFM+FIX+PASSING
                    %每个间距对应的周围车辆速度
                    speed_real = mean(Trj_real(p_index_real(2:end),6));
                    speed_sfm = mean(Trj_sfm(p_index_sfm(2:end),6));
                    speed_sfm_fix_passing = mean(Trj_MSFMcomfortzone(p_index_sfm_fix_passing(2:end),6));
                    
                    L_REAL(c_real:c_real+length(inter_length_real)-1,:) = [inter_length_real speed_real*ones(length(inter_length_real),1)];
                    L_SFM(c_sfm:c_sfm+length(inter_length_sfm)-1,:) = [inter_length_sfm speed_sfm*ones(length(inter_length_sfm),1)];
                    L_SFM_FIX_PASSING(c_sfm_fix_passing:c_sfm_fix_passing+length(inter_length_sfm_fix_passing)-1,:) = [inter_length_sfm_fix_passing speed_sfm_fix_passing*ones(length(inter_length_sfm_fix_passing),1)];
                    c_real = c_real+length(inter_length_real);
                    c_sfm = c_sfm+length(inter_length_sfm);
                    c_sfm_fix_passing = c_sfm_fix_passing+length(inter_length_sfm_fix_passing);
                end
            end
        end
    end
  % 间距直方图分布
    L_REAL = L_REAL (1:length(L_REAL),:);
    L_SFM = L_SFM (1:length(L_SFM),:);
    L_SFM_FIX_PASSING = L_SFM_FIX_PASSING (1:length(L_SFM_FIX_PASSING),:);
    L_REAL = L_REAL(L_REAL (:,1)>0.8,:);
    L_SFM = L_SFM(L_SFM (:,1)>0.8,:);
    L_SFM_FIX_PASSING = L_SFM_FIX_PASSING(L_SFM_FIX_PASSING (:,1)>0.8,:);
    %%数据分组
    % X = [0.5 2 3.5 5 6.5 8 9.5 11 12.5 14];%间距m
    X = [1 2 3 4 5 6 7 8 9 10 ];%间距1m
    numgroup = zeros(11,3);%行为组数
    % 第一列REAL，第二列SFM，第三列SFM_FIX_PASSING
    numgroup(1,1) = sum(L_REAL(:,1)<X(1));
    numgroup(1,2) = sum(L_SFM(:,1)<X(1));
    numgroup(1,3) = sum(L_SFM_FIX_PASSING(:,1)<X(1));
    numgroup(11,1) = sum(L_REAL(:,1)>=X(10));
    numgroup(11,2) = sum(L_SFM(:,1)>=X(10));
    numgroup(11,3) = sum(L_SFM_FIX_PASSING(:,1)>=X(10));
    for i = 2:10
        numgroup(i,1) = sum(L_REAL(:,1)<X(i)&L_REAL(:,1)>=X(i-1));
        numgroup(i,2) = sum(L_SFM(:,1)<X(i)&L_SFM(:,1)>=X(i-1));%
        numgroup(i,3) = sum(L_SFM_FIX_PASSING(:,1)<X(i)&L_SFM_FIX_PASSING(:,1)>=X(i-1));
    end
    for i =1:3
        numgroup(:,i) = numgroup(:,i)./sum(numgroup(:,i));
    end
    figure;
    h = bar(numgroup(1:end,:),'grouped');
%     h = bar(numgroup(1:end,[1 3]),'grouped');
    title 'The distribution of bicycle interspace';
    % set(gca,'XTicklabel',{'0.5' '2' '3.5' '5' '6.5' '8' '9.5' '11' '12.5' '14' '15.5'});
    set(gca,'XTicklabel',{ '<1' '1.5' '2.5' '3.5' '4.5' '5.5' '6.5' '7.5' '8.5' '9.5' '>10' });
    h = legend('FIELD','SFM','MYMODLE');
    set(h, 'Box', 'off');
    set(h,'Interpreter','none')
    xlabel('Interspace among bicycles/m'); ylabel('Frequency');
%     saveas(gcf,'仿真时间帧车辆间间距','fig');
           
    % 间距-速度关系分布   有什么明确的实际意义吗？好像并没有什么意义！！！
    L_REAL = L_REAL (1:length(L_REAL),:);
    L_SFM = L_SFM (1:length(L_SFM),:);
    L_SFM_FIX_PASSING = L_SFM_FIX_PASSING (1:length(L_SFM_FIX_PASSING),:);
    %%数据分组
    X = [0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5];%间距m
    speedgroup = zeros(15,4);%行为组数
    % 第一列REAL，第二列SFM，第三列SFM_FIX，第四列SFM_FIX_PASSING
    speedgroup(1,1) = mean(L_REAL(L_REAL(:,1)<X(1)+0.5,2));
    speedgroup(1,2) = mean(L_SFM(L_SFM(:,1)<X(1)+0.5,2));
    speedgroup(1,3) = mean(L_SFM_FIX_PASSING(L_SFM_FIX_PASSING(:,1)<X(1)+0.5,2));
    speedgroup(15,1) = mean(L_REAL(L_REAL(:,1)>=X(15)+0.5,2));
    speedgroup(15,2) = mean(L_SFM(L_SFM(:,1)>=X(15)+0.5,2));
    speedgroup(15,3) = mean(L_SFM_FIX_PASSING(L_SFM_FIX_PASSING(:,1)>=X(15)+0.5,2));
    for i = 2:14
    speedgroup(i,1) = mean(L_REAL(L_REAL(:,1)<X(i)+0.5&L_REAL(:,2)>=X(i-1)+0.5,2));
    speedgroup(i,2) = mean(L_SFM(L_SFM(:,1)<X(i)+0.5&L_SFM(:,2)>=X(i-1)+0.5,2));
    speedgroup(i,3) = mean(L_SFM_FIX_PASSING(L_SFM_FIX_PASSING(:,1)<X(i)+0.5&L_SFM_FIX_PASSING(:,2)>=X(i-1)+0.5,2));
    end
    figure;
    hold on
    plot(X,speedgroup(:,1),'bo');
    plot(X,speedgroup(:,2),'ks');
    plot(X,speedgroup(:,3),'r^');
    title '车辆速度均值与间距的关系';
    set(gca,'XTicklabel',{'0.5' '1.5' '2.5' '3.5' '4.5' '5.5' '6.5' '7.5' '8.5' '9.5' '10.5' '11.5' '12.5' '13.5' '14.5'});
    h = legend('FIELD','SFM','MSFM');
%     h = legend('FIELD','MYMODEL');    
%     set(h, 'Box', 'off');
%     set(h,'Interpreter','none')
    xlabel('车辆间距/m'); ylabel('平均速度m/s');
    % % saveas(gcf,'仿真时间帧车辆间间距','jpg');
    axis([0.5 13 4 10]);
  grid;
       saveas(gcf,'车辆速度均值与间距的关系','fig');  

% % 超车参数分布   
% figure
% xbins = 0:5:30;
% % histogram(LY_FIELD,xbins);
% h = histogram(LY_FIELD,'Normalization','probability');
% h.BinEdges = [0:5:30];
% hold on
% % histogram(LY_MYMODEL,xbins);
% h = histogram(LY_MYMODEL,'Normalization','probability');
% h.BinEdges = [0:5:30];
% %  saveas(gcf,['超车起始阶段纵向间距分布图','.fig']);
% figure
% xbins = 1:0.5:5;
% % histogram(CLEARENCE_FIELD,xbins);
% h = histogram(CLEARENCE_FIELD,'Normalization','probability');
% h.BinEdges = [1:0.5:5];
% hold on
% % histogram(CLEATENCE_MYMODEL,xbins);
% h = histogram(CLEATENCE_MYMODEL,'Normalization','probability');
% h.BinEdges = [1:0.5:5];
% %   saveas(gcf,['超车起始阶段横向净空分布图','.fig']); 
% figure
% xbins = 0:0.1:3;
% % histogram(LX_FIELD,xbins);
% h = histogram(LX_FIELD,'Normalization','probability');
% h.BinEdges = xbins;
% hold on
% % histogram(LX_MYMODEL,xbins);
% h = histogram(LX_MYMODEL,'Normalization','probability');
% h.BinEdges = xbins;
% %    saveas(gcf,['超车并行横向间距时刻分布图','.fig']);
%    