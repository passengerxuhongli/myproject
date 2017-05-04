%% comfort zone theory
clear;
clc;
bike = xlsread('link-ebike.xlsx','ebike');
%% function 1 : 原始数据时间和格式处理
[start_state] = initdata(bike);
%==============
Tmin=min(start_state(:,2));
Tmax = max(start_state(:,2));
vehicle=start_state(abs(Tmin-start_state(:,2))<0.06,:); 
data =start_state(abs(Tmin-start_state(:,2))<0.06,:); 
tua=3.8737;
A2=0.3627;
B2=34.5654;
lambda = 0.13; 
w=0.4; l=0.8; 
epslow = sqrt(l^2-w^2)/l; 
%% 中间变量和最终变量
information = [];labeltag =[];%舒适空间&行为选择信息记录矩阵
overtaking = [];%记录超车的编号，时刻，纵向间距，横向间距，横向净空
taglastA = zeros(length(vehicle(:,2)),5); 
taglastB = zeros(length(vehicle(:,2)),5);
tag = 0;%超车标签
retag = 0;%超车消失横向速度瞬时间变为0 
for T=Tmin:0.12:Tmax
    %% function 2：discomfort zone-- moving with desired speed
    [a_driving] = drivingforce(vehicle);
    %==========================
     forcebike_bike = zeros(length(vehicle(:,2)),2);
     forceboundary = zeros(length(vehicle(:,2)),2);%边界力 
      %行为标签归零
     TAGA = zeros(length(vehicle(:,2)),2); %左，右各自的标签
     TAGB = zeros(length(vehicle(:,2)),2);
     tap = zeros(length(vehicle(:,2)),7);
     pa_ad = zeros(length(vehicle(:,2)),2);%前车避让横向调整   
       
    for i = 1:size(vehicle,1)
        TAGA(i,:) = [0 0]; 
        TAGB(i,:) = [0 0];
        
        choicetag = zeros(length(vehicle(:,2)),30);
        subject_vehicle = vehicle(vehicle(:,1)==vehicle(i,1),:);
        subject_number = subject_vehicle(1,1);          
        
       %%  function 3 ：舒适空间横纵向边界-------每辆车的舒适空间 
        [comfort_lateral,interdistance,LX,LY,LYB,leadvehicle] = comfortzonefeature(subject_vehicle,vehicle ); 
        
        if exist('leadvehicle','var')&~isempty(leadvehicle)         
            tapleadvehicle(vehicle(i,1),:)=vehicle(vehicle(:,1)==leadvehicle(1,1),:);
        else if exist('tapleadvehicle','var')&size(tapleadvehicle,1)>=i
                try
                    tapleadvehicle(vehicle(i,1),:)=vehicle(vehicle(:,1)==tapleadvehicle(i,1),:);
                end%记录有史以来离我最近的时刻的前车
            end
        end
        overtakeback = [];  %为空表示该车没有前车，如不在自由行使，则处在返回阶段        
         if ~isempty(leadvehicle)
             %记录舒适空间信息
             comfort = [subject_vehicle(1,1),leadvehicle(1,1),T,comfort_lateral,LX,LY,LYB];
            %% function 4 : mid-term model path of choice
    %舒适空间内选择行为模式    
              TLA = taglastA(taglastA(:,1)==subject_vehicle(1,1) &...
                              taglastA(:,2)==leadvehicle(1,1),4:5);
              TLB = taglastB(taglastB(:,1)==subject_vehicle(1,1) &...
                              taglastB(:,2)==leadvehicle(1,1),4:5); 
                overtakeback = 0;%该车有前车，一定不会之在返回阶段
              if isempty(TLA)
                 TLA = [0 0];
              end
              if isempty(TLB)
                 TLB = [0 0];
              end              
             if  sum(TLB,2)==0
                [tagA,tagB,tag_rang,tag_emergency,pathchoice_dataset]=...
                             pathchiocerules(subject_vehicle,vehicle,LX,LY,LYB,leadvehicle,TLA,TLB) ;%L:1 S:2 R:3，不起作用：tag=0                  
             else
               tagB = taptagB(1,vehicle(i,1));%查询中间变量中本车的tagB的值
               tagA =0;
             end
             if tagA == 1
                TAGA(i,:) = [1,0];  
             end
             if tagA == 2 | tagA == 0
                TAGA(i,:) = [0,0]; 
             end   
             if tagA == 3
                TAGA(i,:) = [0,1]; 
             end      
             if tagB == 1
                TAGB(i,:) = [1,0];
             end
             if tagB == 2 | tagB == 0
                TAGB(i,:) = [0,0];
             end   
             if tagB == 3
                TAGB(i,:) = [0,1];
            end   
         else
             comfort = [subject_vehicle(1,1),-1,T,comfort_lateral,LX,LY];        
         end
             information = [information;comfort];     
             tap(i,:) = [comfort(:,1:3),TAGA(i,:),TAGB(i,:)];            
     %% function 5 边界约束函数
      r1 = ((0.5*comfort_lateral)+(1.9*(0.6+interdistance)/0.6)*0.5)/2;%动态的      
      dis_mod_b = vehicle(i,3)-2;%到内边界的距离
      dis_mod1_b = 6-vehicle(i,3);%外部边界
      forceboundary(i,1) = A2*exp(-(dis_mod_b-r1)/B2);
      forceboundary(i,2) = -(A2)*exp(-(dis_mod1_b-r1)/B2);   
     %% function 6 行为执行模型
if ~isempty(leadvehicle)  
     if tagA == 0 & tagB ==0
         % 无行为变化
         forcebike_bike(i,:) = 0;
     end  
     if tagA ==1%向左调整
         %相对速度变化趋势线调整
         [newacceleration] = overtakeapproving(subject_vehicle,leadvehicle);
         forcebike_bike(i,:) = newacceleration;
         a_driving(i,:) = 0;
         forceboundary(i,:) = 0;         
     end
     if tagA ==2%保持直行
         %车辆建无作用，无行为变化
         forcebike_bike(i,:) = 0;
     end
        %========
     if tagA ==3%向右调整
         %相对速度变化趋势线调整
         [newacceleration] = overtakeapproving(subject_vehicle,leadvehicle);
         forcebike_bike(i,:) = newacceleration;
         a_driving(i,:) = 0;
         forceboundary(i,:) = 0;         
     end
        %========  
     if tagB ==1%向左超车
         %相对速度变化趋势线调整
         [newacceleration] = overtakeapproving(subject_vehicle,leadvehicle);
         forcebike_bike(i,:) = newacceleration;
         a_driving(i,:) = 0;
         forceboundary(i,:) = 0;         
     end  
        %========
     if tagB ==2%执行跟车
         forcebike_bike(i,1) = 0.5*(leadvehicle(1,5)-subject_vehicle(1,5));
         forcebike_bike(i,2) = 0.5*(leadvehicle(1,6)-subject_vehicle(1,6))./(leadvehicle(1,4)-subject_vehicle(1,4));%跟车加速度
         a_driving(i,:) = 0;
     end
        %========  
     if tagB ==3%向右超车
         %相对速度变化趋势线调整
         [newacceleration] = overtakeapproving(subject_vehicle,leadvehicle);
         forcebike_bike(i,:) = newacceleration;
         a_driving(i,:) = 0;    
         forceboundary(i,:) = 0;      
     end     
        %========
        if tag_rang==1
    % 判断前车周围是否有条件调整横向位置避让后车
    % ===================
             slvehicle = vehicle(abs(vehicle(:,4)-leadvehicle(1,4))<=5,:);
             slvehicle = slvehicle(slvehicle(:,1)~=leadvehicle(1,1),:);
            if ~isempty(slvehicle)
                %选择可用空间让行
                  pa_ad(i,2)= 0.02;  %具体取正值还是负值应根据可用空间的方向选取   %0.01
                  pa_ad(i,1)= leadvehicle(1,1);  
            else%否则向右调整让行
               pa_ad(i,2)= -0.02;  
               pa_ad(i,1)= leadvehicle(1,1); 
            end       
        end           
     if tag_emergency==1 %紧急情况
         forcebike_bike(i,1) = 0;
         forcebike_bike(i,2) = -0.5;  %0.8
         a_driving(i,:) = 0;    
     end  
          
else  %若前车为空       
        if ~exist('tagB','var')
            tagB = 0;
        end       
      if isempty(overtakeback)& (tagB~=0 & tagB~=2)%无前车的状态只能为自由行使/返回，tagB不为零和2表明该车处于返回阶段
       %the distance between subject_vehicle and leadvehicle      
       distance = subject_vehicle(1,4)-tapleadvehicle(i,4);       
                 if distance>0&distance<5  %返回原有车道的机制
                       newspeed = zeros(1,2);
                       delta_speed_lateral = -0.012*distance+0.23;
                       delta_speed_longit = 0.11*distance+1.5;
                       newspeed(1,1) = tapleadvehicle(i,5)+delta_speed_lateral;
                       newspeed(1,2) = tapleadvehicle(i,6)+delta_speed_longit;
                       newacceleration = (newspeed-subject_vehicle(1,5:6))./0.12;     
                       forcebike_bike(i,:) = newacceleration;
                       a_driving(i,:) = 0;
                       forceboundary(i,:) = 0;  
                 end
      end     
            
end
         taptagB(vehicle(i,1))=tagB;   %给中间变量记录当前车辆这一时刻的tagB,供下一时刻的循环使用   
         if ~isempty(leadvehicle)
            choicetag(i,:) = [tagA,tagB,tag_rang,tag_emergency,pathchoice_dataset,forceboundary(i,:),forcebike_bike(i,:),a_driving(i,:)];
         else
            choicetag(i,:) = [zeros(1,24),forceboundary(i,:),forcebike_bike(i,:),a_driving(i,:)];  
         end  
         
    end

 
pa_ad = pa_ad(pa_ad(:,1)~=0,:);
 
 % 更新
     %更新时间,期望目的地和期望速度不更新。
     vehicle(:,2)= vehicle(:,2)+0.12;
     %更新坐标,t=0.12s
     vehicle(:,3) =  vehicle(:,3) + 1/2*vehicle(:,7)*0.12^2 + vehicle(:,5)*0.12;
     vehicle(:,4) =  vehicle(:,4) + 1/2*vehicle(:,8)*0.12^2 + vehicle(:,6)*0.12;
     %更新速度
     vehicle(:,5) = vehicle(:,5) + vehicle(:,7)*0.12;
     vehicle(:,6) = vehicle(:,6) + vehicle(:,8)*0.12;
     for j=1:length(vehicle(:,2))
if sum(tap(j,6:7),2) ==0%无超车时 
     if abs(vehicle(j,5)) <0.5
          vehicle(j,5) = vehicle(j,5);
     else 
         vehicle(j,5) = (0.5+0.1*rand(1,1))*sign(vehicle(j,5));
     end   
     if vehicle(j,6)<0
        vehicle(j,6) = 0; 
     end
     if vehicle(j,6)>1.1*vehicle(j,11)
        vehicle(j,6) = 1.1*vehicle(j,11); 
     end
    %加速度更新    
    %前车让行加速度更新
%   =============  
    if ~isempty(choicetag(choicetag(:,3)==1,6))&~isempty(pa_ad)
        inter_pa_ad = pa_ad(pa_ad(:,1)==choicetag(choicetag(:,3)==1,6),2);
        if ~isempty(inter_pa_ad)
        vehicle(vehicle(:,1)==choicetag(choicetag(:,3)==1,6),7) = unique( pa_ad(pa_ad(:,1)==choicetag(choicetag(:,3)==1,6),2));   
        end
    end
%     tag_emergency=2;%避让结束标志
%  ==============   
         vehicle(:,7) = forcebike_bike(:,1)+a_driving(:,1)+forceboundary(:,2)+forceboundary(:,1);
         acc8 = forcebike_bike(:,2)+a_driving(:,2);
         acc8(acc8>2) = 1.0;
         acc8(acc8<-2) = -1.5;
         vehicle(:,8) = acc8;
     if abs(vehicle(j,7))<0.6
        vehicle(j,7) = vehicle(j,7);
     else
       vehicle(j,7) = 0.6*sign(vehicle(j,7))+0.1*rand(1,1);%0.8
     end
     if vehicle(j,8)>=0
     if sqrt(vehicle(j,7).^2+vehicle(j,8).^2)<2+0.1*rand(1,1)
         vehicle(j,8) = vehicle(j,8);
     else
          vehicle(j,8) = sqrt(4-vehicle(j,7).^2);
     end
     else if sqrt(vehicle(j,7).^2+vehicle(j,8).^2)<3.5+0.1*rand(1,1)
         vehicle(j,8) = vehicle(j,8);
         else
         vehicle(j,8) = -sqrt(9-vehicle(j,7).^2);
         end
     end
else
    
    if abs(vehicle(j,5)) <0.5%1
          vehicle(j,5) = vehicle(j,5);
     else 
         vehicle(j,5) = (0.5+0.1*rand(1,1))*sign(vehicle(j,5));
     end
     if vehicle(j,6)<0
        vehicle(j,6) = 0; 
     end
     if vehicle(j,6)>1.2*vehicle(j,11)
        vehicle(j,6) = 1.2*vehicle(j,11); 
     end
    %前车让行加速度更新
%   =============  
    if ~isempty(choicetag(choicetag(:,3)==1,6))
    vehicle(vehicle(:,1)==choicetag(choicetag(:,3)==1,6),7) = unique( pa_ad(pa_ad(:,1)==choicetag(choicetag(:,3)==1,6),2));  
    end
%         tag_emergency=2;%避让结束标志
%  ==============   
    %加速度更新
         vehicle(:,7) = forcebike_bike(:,1)+a_driving(:,1)+forceboundary(:,2)+forceboundary(:,1);
         acc8 = forcebike_bike(:,2)+a_driving(:,2);
         acc8(acc8>2) = 1.0;
         acc8(acc8<-2) = -1.5;
         vehicle(:,8) = acc8;
     if abs(vehicle(j,7))<1.24
        vehicle(j,7) = vehicle(j,7);
     else
       vehicle(j,7) = 0.6*sign(vehicle(j,7))+0.1*rand(1,1);
     end
     if vehicle(j,8)>=0
     if sqrt(vehicle(j,7).^2+vehicle(j,8).^2)<2+0.1*rand(1,1)
         vehicle(j,8) = vehicle(j,8);
     else
          vehicle(j,8) = sqrt(4-vehicle(j,7).^2);
     end
     else if sqrt(vehicle(j,7).^2+vehicle(j,8).^2)<3.5+0.1*rand(1,1)
         vehicle(j,8) = vehicle(j,8);
         else
         vehicle(j,8) = -sqrt(9-vehicle(j,7).^2);
         end
     end  
end
     end                  
%========
    % 下一仿真秒，系统车辆更新    
       new_vehicle = start_state(abs(T+0.12-start_state(:,2))<0.06,:);%新生成车辆的编号
       vehicle = [vehicle;new_vehicle];%系统中现有车辆
       data = [data; vehicle];    
 %保留上一秒的选择
if ~isempty(vehicle) & ~isempty(tap) 
A = zeros(size(vehicle,1),5);
B = zeros(size(vehicle,1),5);
taglastA = [A(1:length(tap(:,1)),:)+tap(:,1:5); A(length(tap(:,1))+1:end,:)];
taglastB = [B(1:length(tap(:,1)),:)+tap(:,[1:3,6:7]); B(length(tap(:,1))+1:end,:)];
  %车辆清除  
   index_delete = find(vehicle(:,4)>80);
   taglastA(index_delete,:) = []; 
   taglastB(index_delete,:) = []; 
   vehicle(vehicle(:,4)>80,:)=[];  
else
taglastA = zeros(size(vehicle,1),5);
taglastB = zeros(size(vehicle,1),5);   
end
  % 行为标签记录
  labeltag = [labeltag;choicetag];  %所有车辆信息
 end

overtaking = labeltag(labeltag(:,1)~=0|labeltag(:,2)~=0,:);%超车车辆信息
 %删除重复点
 data_final = [];num = bike(~isnan(bike(:,1)),1);
for n =1:length(num)
data1 = data(data(:,1)==num(n),:);
for m = 1:length(data1(:,2))
    tag = find(data1(:,2)==data1(m,2));
    if length(tag)~=1
     data1(tag(2:end),:)=nan;
    end
end
data1(isnan(data1(:,1)),:) = [];
data_final = [data_final;data1];
end
% xlswrite('MSFMcomfortzone_评价基础数据',data_final);
% xlswrite('MSFMcomfortzoneinformation_finalmlogit1',information);
% xlswrite('MSFMcomfortzonedensity-overtaking_finalmlogit1',overtaking);
xlswrite('MSFMcomfortzone_mymodel',data_final);
%% 画图
index = find(~isnan(bike(:,1)));
figure 

for i=1:size(index,1)-1
plot(bike(index(i)+1:index(i+1)-1,4),bike(index(i)+1:index(i+1)-1,3),'b');
hold on;
plot( data_final(data_final(:,1)==bike(index(i),1),4),data_final(data_final(:,1)==bike(index(i),1),3),'r'); 
% pause(0.5);
% hold off;
end
axis([-10 100 0 8]);
grid;

    