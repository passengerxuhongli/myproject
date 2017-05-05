%% 仿真可视化
%读取数据
clear
clc
Trj_real = xlsread('评价基础数据.xlsx',1);
Trj_MSFM = xlsread('评价基础数据.xlsx',4);
data_final = Trj_real;
data_final = data_final(data_final(:,4)<=75,:);
Tmin = min(data_final(:,2));
Tmax = max(data_final(:,2));
aviobj = VideoWriter('C:\Users\Administrator\Desktop\comfortzoneofbicycle\evaluation\版本二\field.avi');
aviobj.FrameRate = 10;
open(aviobj);
scrsz = get(0,'ScreenSize');
for T=Tmin:0.12:Tmax%672:0.12:980
% 初始化
% 路网加载
hh=figure('numbertitle','off','name',...
'车辆运动图','color',([0.15 0.15 0.15]),'Position',[1 1 scrsz(3) scrsz(4)]);  %[0 0 0]是黑色
hold on;
axis off;
      text(25,10,'电动车超越自行车车辆运动可视化',...
'fontname','楷体','fontsize',18,'FontWeight','demi','Color','w');    
% 车辆加载
      current_data = data_final(abs(data_final(:,2)-T)<0.06,:);    
 %========防止出现两个一样的车辆编号======     
for i=1:length(current_data(:,1))-1
 for j=i+1:length(current_data(:,1))
  if current_data(j,1)==current_data(i,1)
   current_data(j,:)=0;%或者其他标识符号。
  end
 end
end
idx=find(current_data(:,1)==0);%与前面的标识符号一致。
current_data(idx,:)=[];%删除标识项
 %========防止出现两个一样的车辆编号======
 
for i= 1:length(current_data(:,1))
theta=0:pi/20:2*pi;%画椭圆
fill(current_data(current_data(:,1)==current_data(i,1),4)+0.9*cos(theta),current_data(current_data(:,1)==current_data(i,1),3)+0.3*sin(theta),'b');%颜色填充
hold on
%% 超车车辆舒适空间椭圆
subject_vehicle = current_data(i,:);
vehicle = current_data;
Q = sum(vehicle(:,4)-subject_vehicle(1,4)>=-5 & vehicle(:,4)-subject_vehicle(1,4)<=5);
density = Q/(4*10);
interdistance = min(0.156*density^(-0.18746),1.6);
comfort_lateral = min(0.6+interdistance,1.6);
lx = 0.8-0.5*subject_vehicle(1,11)*sin(pi-2*acos(9.8*1.0*tan(pi*15/180)/(subject_vehicle(1,11).^2)));
LX = lx+comfort_lateral;
%前车
vehicle1 =  vehicle(vehicle(:,1)~=subject_vehicle(1,1),:);
leadvehicle = vehicle1(abs(vehicle1(:,3)-subject_vehicle(1,3))<=LX & vehicle1(:,4)-subject_vehicle(1,4)>=0,:);
    v_subjectvehicle = subject_vehicle(1,6);
if ~isempty(leadvehicle)
    leadvehicle = leadvehicle(leadvehicle(:,4)==min(leadvehicle(:,4)),:);
    Q1 = sum(abs(vehicle(:,4)-leadvehicle(1,4))<=2.5);
    k_leadvehicle = Q1/(4*5);
end
if  ~isempty(leadvehicle) 
     v_leadvehicle = leadvehicle(1,6);
  if v_subjectvehicle>=v_leadvehicle
     LY = abs(15.3911+0.5*v_subjectvehicle+0.1784*v_subjectvehicle.^2-0.4234*v_leadvehicle.^2-56.1606*k_leadvehicle);
  else 
     LY =  v_subjectvehicle*0.5+v_subjectvehicle^2/(2*2);   
  end
  
%判断超车，若超车则前车（被超车）填充变为红色，本车变为绿色。并显示超车对编号
if leadvehicle(1,4)-subject_vehicle(1,4)<=LY & abs(leadvehicle(1,3)-subject_vehicle(1,3))<=LX
%本车编号
 text(subject_vehicle(1,4),subject_vehicle(1,3),num2str(subject_vehicle(1,1)),...
'fontname','楷体','fontsize',12,'FontWeight','demi','Color','r');    
%前车编号
 text(leadvehicle(1,4),leadvehicle(1,3),num2str(leadvehicle(1,1)),...
'fontname','楷体','fontsize',12,'FontWeight','demi','Color','g'); 
%变更颜色
fill(subject_vehicle(1,4)+0.9*cos(theta),subject_vehicle(1,3)+0.3*sin(theta),'r');%颜色填充
fill(leadvehicle(1,4)+0.9*cos(theta),leadvehicle(1,3)+0.3*sin(theta),'g');%颜色填充
end  
   
else 
LY =  v_subjectvehicle*0.5+v_subjectvehicle^2/(2*2);   
end
%舒适空间展示（动态的）
theta=-0.5*pi:pi/20:0.5*pi;%画椭圆
plot(-1+subject_vehicle(1,4)+LY*cos(theta), subject_vehicle(1,3)+LX*sin(theta),'g-','linewidth',1);%自行车0.9 0.3  1.8 0.6
hold on 
theta=0.5*pi:pi/20:1.5*pi;%画椭圆
plot(-1+subject_vehicle(1,4)+LX*cos(theta), subject_vehicle(1,3)+LX*sin(theta),'g-','linewidth',1);%自行车0.9 0.3  1.8 0.6
% axis equal

end
axis equal;
linestart=-15;
lineend=110;%路段长度
plot(linestart:1:lineend,6+zeros(1,(lineend-linestart)+1),'w','linewidth',2);
hold on
plot(linestart:1:lineend,2+zeros(1,(lineend-linestart)+1),'w','linewidth',2);
plot(linestart:1:lineend,4+zeros(1,(lineend-linestart)+1),'w-.','linewidth',.2);
plot(linestart:1:lineend,3+zeros(1,(lineend-linestart)+1),'w-.','linewidth',0.2);
plot(linestart:1:lineend,5+zeros(1,(lineend-linestart)+1),'w-.','linewidth',.2);
    frame = getframe(gcf);
    writeVideo(aviobj,frame);
    delete(gcf);
end 
 close(aviobj); 

%% ========== MSFM ===============
data_final = Trj_MSFM;
data_final = data_final(data_final(:,4)<=75,:);
Tmin = min(data_final(:,2));
Tmax = max(data_final(:,2));
aviobj = VideoWriter('C:\Users\Administrator\Desktop\comfortzoneofbicycle\evaluation\版本二\mymodel.avi');
aviobj.FrameRate = 10;
open(aviobj);
scrsz = get(0,'ScreenSize');
for T=Tmin:0.12:Tmax 
% 初始化
% 路网加载
hh=figure('numbertitle','off','name',...
'车辆运动图','color',([0.15 0.15 0.15]),'Position',[1 1 scrsz(3) scrsz(4)]);  %[0 0 0]是黑色
hold on;
axis off;
      text(25,10,'电动车超越自行车车辆运动可视化',...
'fontname','楷体','fontsize',18,'FontWeight','demi','Color','w');    
% 车辆加载
      current_data = data_final(abs(data_final(:,2)-T)<0.06,:);    
 %========防止出现两个一样的车辆编号======     
for i=1:length(current_data(:,1))-1
 for j=i+1:length(current_data(:,1))
  if current_data(j,1)==current_data(i,1)
   current_data(j,:)=0;%或者其他标识符号。
  end
 end
end
idx=find(current_data(:,1)==0);%与前面的标识符号一致。
current_data(idx,:)=[];%删除标识项
 %========防止出现两个一样的车辆编号======
 
for i= 1:length(current_data(:,1))
theta=0:pi/20:2*pi;%画椭圆
fill(current_data(current_data(:,1)==current_data(i,1),4)+0.9*cos(theta),current_data(current_data(:,1)==current_data(i,1),3)+0.3*sin(theta),'b');%颜色填充
hold on
%% 超车车辆舒适空间椭圆
subject_vehicle = current_data(i,:);
vehicle = current_data;
Q = sum(vehicle(:,4)-subject_vehicle(1,4)>=-5 & vehicle(:,4)-subject_vehicle(1,4)<=5);
density = Q/(4*10);
interdistance = min(0.156*density^(-0.18746),1.6);
comfort_lateral = min(0.6+interdistance,1.6);
lx = 0.8-0.5*subject_vehicle(1,11)*sin(pi-2*acos(9.8*1.0*tan(pi*15/180)/(subject_vehicle(1,11).^2)));
LX = lx+comfort_lateral;
%前车
vehicle1 =  vehicle(vehicle(:,1)~=subject_vehicle(1,1),:);
leadvehicle = vehicle1(abs(vehicle1(:,3)-subject_vehicle(1,3))<=LX & vehicle1(:,4)-subject_vehicle(1,4)>=0,:);
    v_subjectvehicle = subject_vehicle(1,6);
if ~isempty(leadvehicle)
    leadvehicle = leadvehicle(leadvehicle(:,4)==min(leadvehicle(:,4)),:);
    Q1 = sum(abs(vehicle(:,4)-leadvehicle(1,4))<=2.5);
    k_leadvehicle = Q1/(4*5);
end
if  ~isempty(leadvehicle) 
     v_leadvehicle = leadvehicle(1,6);
  if v_subjectvehicle>=v_leadvehicle
     LY = abs(15.3911+0.5*v_subjectvehicle+0.1784*v_subjectvehicle.^2-0.4234*v_leadvehicle.^2-56.1606*k_leadvehicle);
  else 
     LY =  v_subjectvehicle*0.5+v_subjectvehicle^2/(2*2);   
  end
  
%判断超车，若超车则前车（被超车）填充变为红色，本车变为绿色。并显示超车对编号
if leadvehicle(1,4)-subject_vehicle(1,4)<=LY & abs(leadvehicle(1,3)-subject_vehicle(1,3))<=LX
%本车编号
 text(subject_vehicle(1,4),subject_vehicle(1,3),num2str(subject_vehicle(1,1)),...
'fontname','楷体','fontsize',12,'FontWeight','demi','Color','r');    
%前车编号
 text(leadvehicle(1,4),leadvehicle(1,3),num2str(leadvehicle(1,1)),...
'fontname','楷体','fontsize',12,'FontWeight','demi','Color','g'); 
%变更颜色
fill(subject_vehicle(1,4)+0.9*cos(theta),subject_vehicle(1,3)+0.3*sin(theta),'r');%颜色填充
fill(leadvehicle(1,4)+0.9*cos(theta),leadvehicle(1,3)+0.3*sin(theta),'g');%颜色填充
end  
   
else 
LY =  v_subjectvehicle*0.5+v_subjectvehicle^2/(2*2);   
end
%舒适空间展示（动态的）
theta=-0.5*pi:pi/20:0.5*pi;%画椭圆
plot(-1+subject_vehicle(1,4)+LY*cos(theta), subject_vehicle(1,3)+LX*sin(theta),'g-','linewidth',1);%自行车0.9 0.3  1.8 0.6
hold on 
theta=0.5*pi:pi/20:1.5*pi;%画椭圆
plot(-1+subject_vehicle(1,4)+LX*cos(theta), subject_vehicle(1,3)+LX*sin(theta),'g-','linewidth',1);%自行车0.9 0.3  1.8 0.6


end
axis equal;
linestart=-15;
lineend=110;%路段长度
plot(linestart:1:lineend,6+zeros(1,(lineend-linestart)+1),'w','linewidth',2);
hold on
plot(linestart:1:lineend,2+zeros(1,(lineend-linestart)+1),'w','linewidth',2);
plot(linestart:1:lineend,4+zeros(1,(lineend-linestart)+1),'w-.','linewidth',.2);
plot(linestart:1:lineend,3+zeros(1,(lineend-linestart)+1),'w-.','linewidth',0.2);
plot(linestart:1:lineend,5+zeros(1,(lineend-linestart)+1),'w-.','linewidth',.2);
    frame = getframe(gcf);
    writeVideo(aviobj,frame);
    delete(gcf);
end 
 close(aviobj); 
