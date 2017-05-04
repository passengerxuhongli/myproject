%% ������ӻ�
%��ȡ����
clear
clc
Trj_real = xlsread('���ۻ�������.xlsx',1);
Trj_MSFM = xlsread('���ۻ�������.xlsx',4);
data_final = Trj_real;
data_final = data_final(data_final(:,4)<=75,:);
Tmin = min(data_final(:,2));
Tmax = max(data_final(:,2));
aviobj = VideoWriter('C:\Users\Administrator\Desktop\comfortzoneofbicycle\evaluation\�汾��\field.avi');
aviobj.FrameRate = 10;
open(aviobj);
scrsz = get(0,'ScreenSize');
for T=Tmin:0.12:Tmax%672:0.12:980
% ��ʼ��
% ·������
hh=figure('numbertitle','off','name',...
'�����˶�ͼ','color',([0.15 0.15 0.15]),'Position',[1 1 scrsz(3) scrsz(4)]);  %[0 0 0]�Ǻ�ɫ
hold on;
axis off;
      text(25,10,'�綯����Խ���г������˶����ӻ�',...
'fontname','����','fontsize',18,'FontWeight','demi','Color','w');    
% ��������
      current_data = data_final(abs(data_final(:,2)-T)<0.06,:);    
 %========��ֹ��������һ���ĳ������======     
for i=1:length(current_data(:,1))-1
 for j=i+1:length(current_data(:,1))
  if current_data(j,1)==current_data(i,1)
   current_data(j,:)=0;%����������ʶ���š�
  end
 end
end
idx=find(current_data(:,1)==0);%��ǰ��ı�ʶ����һ�¡�
current_data(idx,:)=[];%ɾ����ʶ��
 %========��ֹ��������һ���ĳ������======
 
for i= 1:length(current_data(:,1))
theta=0:pi/20:2*pi;%����Բ
fill(current_data(current_data(:,1)==current_data(i,1),4)+0.9*cos(theta),current_data(current_data(:,1)==current_data(i,1),3)+0.3*sin(theta),'b');%��ɫ���
hold on
%% �����������ʿռ���Բ
subject_vehicle = current_data(i,:);
vehicle = current_data;
Q = sum(vehicle(:,4)-subject_vehicle(1,4)>=-5 & vehicle(:,4)-subject_vehicle(1,4)<=5);
density = Q/(4*10);
interdistance = min(0.156*density^(-0.18746),1.6);
comfort_lateral = min(0.6+interdistance,1.6);
lx = 0.8-0.5*subject_vehicle(1,11)*sin(pi-2*acos(9.8*1.0*tan(pi*15/180)/(subject_vehicle(1,11).^2)));
LX = lx+comfort_lateral;
%ǰ��
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
  
%�жϳ�������������ǰ����������������Ϊ��ɫ��������Ϊ��ɫ������ʾ�����Ա��
if leadvehicle(1,4)-subject_vehicle(1,4)<=LY & abs(leadvehicle(1,3)-subject_vehicle(1,3))<=LX
%�������
 text(subject_vehicle(1,4),subject_vehicle(1,3),num2str(subject_vehicle(1,1)),...
'fontname','����','fontsize',12,'FontWeight','demi','Color','r');    
%ǰ�����
 text(leadvehicle(1,4),leadvehicle(1,3),num2str(leadvehicle(1,1)),...
'fontname','����','fontsize',12,'FontWeight','demi','Color','g'); 
%�����ɫ
fill(subject_vehicle(1,4)+0.9*cos(theta),subject_vehicle(1,3)+0.3*sin(theta),'r');%��ɫ���
fill(leadvehicle(1,4)+0.9*cos(theta),leadvehicle(1,3)+0.3*sin(theta),'g');%��ɫ���
end  
   
else 
LY =  v_subjectvehicle*0.5+v_subjectvehicle^2/(2*2);   
end
%���ʿռ�չʾ����̬�ģ�
theta=-0.5*pi:pi/20:0.5*pi;%����Բ
plot(-1+subject_vehicle(1,4)+LY*cos(theta), subject_vehicle(1,3)+LX*sin(theta),'g-','linewidth',1);%���г�0.9 0.3  1.8 0.6
hold on 
theta=0.5*pi:pi/20:1.5*pi;%����Բ
plot(-1+subject_vehicle(1,4)+LX*cos(theta), subject_vehicle(1,3)+LX*sin(theta),'g-','linewidth',1);%���г�0.9 0.3  1.8 0.6
% axis equal

end
axis equal;
linestart=-15;
lineend=110;%·�γ���
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
aviobj = VideoWriter('C:\Users\Administrator\Desktop\comfortzoneofbicycle\evaluation\�汾��\mymodel.avi');
aviobj.FrameRate = 10;
open(aviobj);
scrsz = get(0,'ScreenSize');
for T=Tmin:0.12:Tmax 
% ��ʼ��
% ·������
hh=figure('numbertitle','off','name',...
'�����˶�ͼ','color',([0.15 0.15 0.15]),'Position',[1 1 scrsz(3) scrsz(4)]);  %[0 0 0]�Ǻ�ɫ
hold on;
axis off;
      text(25,10,'�綯����Խ���г������˶����ӻ�',...
'fontname','����','fontsize',18,'FontWeight','demi','Color','w');    
% ��������
      current_data = data_final(abs(data_final(:,2)-T)<0.06,:);    
 %========��ֹ��������һ���ĳ������======     
for i=1:length(current_data(:,1))-1
 for j=i+1:length(current_data(:,1))
  if current_data(j,1)==current_data(i,1)
   current_data(j,:)=0;%����������ʶ���š�
  end
 end
end
idx=find(current_data(:,1)==0);%��ǰ��ı�ʶ����һ�¡�
current_data(idx,:)=[];%ɾ����ʶ��
 %========��ֹ��������һ���ĳ������======
 
for i= 1:length(current_data(:,1))
theta=0:pi/20:2*pi;%����Բ
fill(current_data(current_data(:,1)==current_data(i,1),4)+0.9*cos(theta),current_data(current_data(:,1)==current_data(i,1),3)+0.3*sin(theta),'b');%��ɫ���
hold on
%% �����������ʿռ���Բ
subject_vehicle = current_data(i,:);
vehicle = current_data;
Q = sum(vehicle(:,4)-subject_vehicle(1,4)>=-5 & vehicle(:,4)-subject_vehicle(1,4)<=5);
density = Q/(4*10);
interdistance = min(0.156*density^(-0.18746),1.6);
comfort_lateral = min(0.6+interdistance,1.6);
lx = 0.8-0.5*subject_vehicle(1,11)*sin(pi-2*acos(9.8*1.0*tan(pi*15/180)/(subject_vehicle(1,11).^2)));
LX = lx+comfort_lateral;
%ǰ��
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
  
%�жϳ�������������ǰ����������������Ϊ��ɫ��������Ϊ��ɫ������ʾ�����Ա��
if leadvehicle(1,4)-subject_vehicle(1,4)<=LY & abs(leadvehicle(1,3)-subject_vehicle(1,3))<=LX
%�������
 text(subject_vehicle(1,4),subject_vehicle(1,3),num2str(subject_vehicle(1,1)),...
'fontname','����','fontsize',12,'FontWeight','demi','Color','r');    
%ǰ�����
 text(leadvehicle(1,4),leadvehicle(1,3),num2str(leadvehicle(1,1)),...
'fontname','����','fontsize',12,'FontWeight','demi','Color','g'); 
%�����ɫ
fill(subject_vehicle(1,4)+0.9*cos(theta),subject_vehicle(1,3)+0.3*sin(theta),'r');%��ɫ���
fill(leadvehicle(1,4)+0.9*cos(theta),leadvehicle(1,3)+0.3*sin(theta),'g');%��ɫ���
end  
   
else 
LY =  v_subjectvehicle*0.5+v_subjectvehicle^2/(2*2);   
end
%���ʿռ�չʾ����̬�ģ�
theta=-0.5*pi:pi/20:0.5*pi;%����Բ
plot(-1+subject_vehicle(1,4)+LY*cos(theta), subject_vehicle(1,3)+LX*sin(theta),'g-','linewidth',1);%���г�0.9 0.3  1.8 0.6
hold on 
theta=0.5*pi:pi/20:1.5*pi;%����Բ
plot(-1+subject_vehicle(1,4)+LX*cos(theta), subject_vehicle(1,3)+LX*sin(theta),'g-','linewidth',1);%���г�0.9 0.3  1.8 0.6


end
axis equal;
linestart=-15;
lineend=110;%·�γ���
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
