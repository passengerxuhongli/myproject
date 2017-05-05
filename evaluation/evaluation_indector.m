
clear
clc
%%
% Trj_real = xlsread('link-ebike.xlsx',1);
% for i= 1:length(index)-1
%    Trj_real(find(Trj_real(:,1)==Trj_real(index(i),1)):find(Trj_real(:,1)==Trj_real(index(i+1),1))-1,1) = Trj_real(index(i),1);
% end%�ѹ켣��Ų�ȫ
% Trj_real(isnan(Trj_real(:,2)),:)=[];%�޳�������
% t0 = min(Trj_real(:,2));%�����ʼʱ��
% Trj_real(:,2) = (Trj_real(:,2)-t0)*24*3600;%ͬһʱ��̶�Ϊ0.12��ı���
% % ��link-ebike.xlsx���ݸ�ʽ��ɡ����ۻ������ݡ��ĸ�ʽ
%��ȡ����
Trj_real = xlsread('���ۻ�������.xlsx',1);
Trj_sfm = xlsread('���ۻ�������.xlsx',2);
Trj_MSFMcomfortzone = xlsread('���ۻ�������.xlsx',3);
index = find(~isnan(Trj_real(:,1)));%real
%% ================ ����Ч������ ==============================
%% ����������������������Ȳ�����ȡ 
%field data ������Ϣ
overtake_number=0;bike_number=[];left_field = 0;times = 0;%��೬�� 
time_overtake = []; longitudinal_distance =[]; lateral_space =[];lateral_distance =[];tip =0 ;%������¼��ǩ
for i=1:max(Trj_real(:,1))%�󳵳���ѭ��
    for T= min(Trj_real(Trj_real(:,1)==i,2)):0.12:max(Trj_real(Trj_real(:,1)==i,2))-0.12%����������ʱ�����
        front_bike = Trj_real(abs(Trj_real(:,2)-T)<0.06,:);%Ѱ������ǰ��         
        back_bike = Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-T)<0.06,:);%Ѱ�ҵ�ǰ��
        front_bike(front_bike(:,1)==i,: )=[]; %�޳�������        
        front_bike1 = Trj_real(abs(Trj_real(:,2)-T-0.12)<0.06,:);%Ѱ����һ��������ǰ��
        back_bike1 = Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-T-0.12)<0.06,:);%Ѱ����һ�����ĺ�
        front_bike1(front_bike1(:,1)==i,: )=[]; %�޳�������
        if ~isempty(front_bike(:,4))
            for j=1:length(front_bike(:,1))
             if  ~isempty( find(front_bike1(:,1)==front_bike(j,1)))
                if  (back_bike(1,4)-front_bike(j,4))<0 &(back_bike1(1,4)-front_bike1(front_bike1(:,1)==front_bike(j,1),4))>=0  %������  %�жϸó����Ƿ�Ϊǰ��&�Ƿ����󳵳�Խǰ������
                    for t=max(min(Trj_real(Trj_real(:,1)==i,2)),min(Trj_real(Trj_real(:,1)==front_bike(j,1),2))):0.12:T
                        zone = Trj_real(Trj_real(:,4)>=Trj_real(Trj_real(:,1)==front_bike(j,1)&abs(Trj_real(:,2)-t)<0.06,4)-2.5&Trj_real(:,4)<=Trj_real(Trj_real(:,1)==front_bike(j,1)&abs(Trj_real(:,2)-t)<0.06,4)+2.5,:);%ͳ������Ϊ5m
                        density = length(zone(abs(zone(:,2)-t)<0.06,1))/(4*5);
                        r1 = abs( 15.3911+0.5*Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,6)+0.1784*Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,6).^2-0.4234*Trj_real(Trj_real(:,1)==front_bike(j,1)&abs(Trj_real(:,2)-t)<0.06,6).^2-56.1606*density ); %����
                        zone1 = Trj_real(Trj_real(:,4)>=Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,4)-5&Trj_real(:,4)<=Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,4)+5,:);%ͳ������Ϊ10m
                        density1 = length(zone1(abs(zone1(:,2)-t)<0.06,1))/(4*10);
                        interdistance = min(0.156*density1^(-0.18746),1.6);
                        comfort_lateral = min(0.8+interdistance,1.6);
                        r2 = 0.8-0.5*Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,6)*sin(pi-2*acos(9.*1.0*tan(pi*15/180)/(Trj_real(Trj_real(:,1)==i&abs(Trj_real(:,2)-t)<0.06,6).^2)))+comfort_lateral;
                        a1=Trj_real(abs(Trj_real(:,2)-t)<0.06&Trj_real(:,1)==i,4);
                        a2=Trj_real(abs(Trj_real(:,2)-t)<0.06&Trj_real(:,1)==front_bike(j,1),4);
                        if abs(a1(1)-a2(1))<=r1
                            a3=Trj_real(abs(Trj_real(:,2)-t)<0.06&Trj_real(:,1)==i,3);
                            a4=Trj_real(abs(Trj_real(:,2)-t)<0.06&Trj_real(:,1)==front_bike(j,1),3);
                            if abs(a3(1)-a4(1))<=r2%�ַ����ʿռ䣿                               
%                                delta_x = abs(front_bike(j,3)-back_bike(1,3));
%                                delta_y = abs(front_bike(j,4)-back_bike(1,4));
%                                 yy = delta_y.^2/r1 + delta_x.^2/r2; 
%                                 if yy<=1 & (front_bike(j,4)-back_bike(1,4)>=0)
                                overtake_number = overtake_number+1;%���������ۼ�
                                time_overtake =[time_overtake; t];%����ʱ��
                                longitudinal_distance =[ longitudinal_distance; abs(a2(1)-a1(1))];%������m
                                lateral_space =[ lateral_space; max(a4(1)-2,6-a4(1))];%������ʼʱǰ����߽�֮��ļ�����ֵ
                                lateral_distance =[lateral_distance; abs( front_bike1(j,3)-back_bike1(1,3))];%�������е㴦��������ࡣ
                                if front_bike(j,3)-back_bike(1,3)>=0
                                    left_field = left_field+1;
                                end
                                tip = 1;
                                bike_number=[bike_number;i];%��������-��
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
    %�������
    times
    overtake_number
    left_field
    xlswrite('��������1',[bike_number,time_overtake,longitudinal_distance,lateral_space,lateral_distance],1);

    %sfm ���������ȳ�����Ϣ
    overtake_number=0;bike_number=[];left_sfm = 0;times = 0;%��೬��
    time_overtake = []; longitudinal_distance =[]; lateral_space =[];lateral_distance =[];tip =0 ;%������¼��ǩ
    for i=1:max(Trj_sfm(:,1))%�󳵳���ѭ��
        for T= min(Trj_sfm(Trj_sfm(:,1)==i,2)):0.12:max(Trj_sfm(Trj_sfm(:,1)==i,2))-0.12%����������ʱ�����
            front_bike = Trj_sfm(abs(Trj_sfm(:,2)-T)<0.06,:);%Ѱ������ǰ��
            back_bike = Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-T)<0.06,:);%Ѱ�ҵ�ǰ��
            front_bike(front_bike(:,1)==i,: )=[]; %�޳�������
            front_bike1 = Trj_sfm(abs(Trj_sfm(:,2)-T-0.12)<0.06,:);%Ѱ����һ��������ǰ��
            back_bike1 = Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-T-0.12)<0.06,:);%Ѱ����һ�����ĺ�
            front_bike1(front_bike1(:,1)==i,: )=[]; %�޳�������
            if ~isempty(front_bike(:,4))
                for j=1:length(front_bike(:,1))
                 if  ~isempty( find(front_bike1(:,1)==front_bike(j,1)))
                    if   (back_bike(1,4)-front_bike(j,4))<0 & (back_bike1(1,4)-front_bike1(front_bike1(:,1)==front_bike(j,1),4))>=0  %������
                        for t= max(min(Trj_sfm(Trj_sfm(:,1)==i,2)),min(Trj_sfm(Trj_sfm(:,1)==front_bike(j,1),2))):0.12:T                                                       
                        zone = Trj_sfm(Trj_sfm(:,4)>=Trj_sfm(Trj_sfm(:,1)==front_bike(j,1)&abs(Trj_sfm(:,2)-t)<0.06,4)-2.5&Trj_sfm(:,4)<=Trj_sfm(Trj_sfm(:,1)==front_bike(j,1)&abs(Trj_sfm(:,2)-t)<0.06,4)+2.5,:);%ͳ������Ϊ5m
                        density = length(zone(abs(zone(:,2)-t)<0.06,1))/(4*5);
                        r1 = abs( 15.3911+0.5*Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,6)+0.1784*Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,6).^2-0.4234*Trj_sfm(Trj_sfm(:,1)==front_bike(j,1)&abs(Trj_sfm(:,2)-t)<0.06,6).^2-56.1606*density ); %����
                        zone1 = Trj_sfm(Trj_sfm(:,4)>=Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,4)-5&Trj_sfm(:,4)<=Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,4)+5,:);%ͳ������Ϊ10m
                        density1 = length(zone1(abs(zone1(:,2)-t)<0.06,1))/(4*10);
                        interdistance = min(0.156*density1^(-0.18746),1.6);
                        comfort_lateral = min(0.6+interdistance,1.6);
                        r2 = 0.8-0.5*Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,6)*sin(pi-2*acos(9.*1.0*tan(pi*15/180)/(Trj_sfm(Trj_sfm(:,1)==i&abs(Trj_sfm(:,2)-t)<0.06,6).^2)))+comfort_lateral;                            
                            a1=Trj_sfm(abs(Trj_sfm(:,2)-t)<0.06&Trj_sfm(:,1)==i,4);
                            a2=Trj_sfm(abs(Trj_sfm(:,2)-t)<0.06&Trj_sfm(:,1)==front_bike(j,1),4);
                            if abs(a1(1)-a2(1))<=r1
                                a3=Trj_sfm(abs(Trj_sfm(:,2)-t)<0.06&Trj_sfm(:,1)==i,3);
                                a4=Trj_sfm(abs(Trj_sfm(:,2)-t)<0.06&Trj_sfm(:,1)==front_bike(j,1),3);
                                if abs(a3(1)-a4(1))<=r2%�ַ����ʿռ䣿
                                    
                                    delta_x = abs(front_bike(j,3)-back_bike(1,3));
                                    delta_y = abs(front_bike(j,4)-back_bike(1,4));
                                    yy = delta_y.^2/r1 + delta_x.^2/r2; 
                                    if yy<=1 & (front_bike(j,4)-back_bike(1,4)>=0)                                    
                                    overtake_number = overtake_number+1;%���������ۼ�
                                    time_overtake =[time_overtake; t];%����ʱ��
                                    longitudinal_distance =[ longitudinal_distance; abs(a2(1)-a1(1))];%������m
                                    lateral_space =[ lateral_space; max(a4(1)-2,6-a4(1))];%������ʼʱǰ����߽�֮��ļ�����ֵ
                                    lateral_distance =[lateral_distance; abs( front_bike1(j,3)-back_bike1(1,3))];%�������е㴦��������ࡣ
                                    if front_bike(j,3)-back_bike(1,3)>=0
                                        left_sfm = left_sfm+1;
                                    end
                                    tip = 1;
                                    bike_number=[bike_number;i];%��������-��
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
    %�������
    times
    overtake_number
    left_sfm
    xlswrite('��������1',[bike_number,time_overtake,longitudinal_distance,lateral_space,lateral_distance],2);

    %msfmcomfortzone ���������ȳ�����Ϣ
    overtake_number=0;bike_number=[];left_msfm = 0;times = 0;%��೬��
    time_overtake = []; longitudinal_distance =[]; lateral_space =[];lateral_distance =[];tip =0 ;%������¼��ǩ
    for i=1:max(Trj_MSFMcomfortzone(:,1))%�󳵳���ѭ��
        for T= min(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,2)):0.12:max(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,2))-0.12%����������ʱ�����
            
            front_bike = Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-T)<0.06,:);%Ѱ������ǰ��
            back_bike = Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-T)<0.06,:);%Ѱ�ҵ�ǰ��
            front_bike(front_bike(:,1)==i,: )=[]; %�޳�������
            front_bike1 = Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-T-0.12)<0.06,:);%Ѱ����һ��������ǰ��
            back_bike1 = Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-T-0.12)<0.06,:);%Ѱ����һ�����ĺ�
            front_bike1(front_bike1(:,1)==i,: )=[]; %�޳�������

            if ~isempty(front_bike(:,4))
                for j=1:length(front_bike(:,1))
                  if  ~isempty( find(front_bike1(:,1)==front_bike(j,1)))
                    if   (back_bike(1,4)-front_bike(j,4))<0 & (back_bike1(1,4)-front_bike1(front_bike1(:,1)==front_bike(j,1),4))>=0  %������
                        for t= max(min(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i,2)),min(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==front_bike(j,1),2))):0.12:T   
                        zone = Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,4)>=Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==front_bike(j,1)...
                            &abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,4)-2.5&Trj_MSFMcomfortzone(:,4)<=Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==front_bike(j,1)...
                            &abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,4)+2.5,:);%ͳ������Ϊ5m
                        density = length(zone(abs(zone(:,2)-t)<0.06,1))/(4*5);
                        r1 = abs( 15.3911+0.5*Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6)+0.1784*Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6).^2-0.4234*Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==front_bike(j,1)&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6).^2-56.1606*density ); %����
                        zone1 = Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,4)>=Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,4)-5&Trj_MSFMcomfortzone(:,4)<=Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,4)+5,:);%ͳ������Ϊ10m
                        density1 = length(zone1(abs(zone1(:,2)-t)<0.06,1))/(4*10);
                        interdistance = min(0.156*density1^(-0.18746),1.6);
                        comfort_lateral = min(0.6+interdistance,1.6);
                        r2 = 0.8-0.5*Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6)*sin(pi-2*acos(9.*1.0*tan(pi*15/180)/(Trj_MSFMcomfortzone(Trj_MSFMcomfortzone(:,1)==i&abs(Trj_MSFMcomfortzone(:,2)-t)<0.06,6).^2)))+comfort_lateral;                            
                            a1=Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-t)<0.06&Trj_MSFMcomfortzone(:,1)==i,4);
                            a2=Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-t)<0.06&Trj_MSFMcomfortzone(:,1)==front_bike(j,1),4);
                            if abs(a1(1)-a2(1))<=r1
                                a3=Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-t)<0.06&Trj_MSFMcomfortzone(:,1)==i,3);
                                a4=Trj_MSFMcomfortzone(abs(Trj_MSFMcomfortzone(:,2)-t)<0.06&Trj_MSFMcomfortzone(:,1)==front_bike(j,1),3);
                                if abs(a3(1)-a4(1))<=r2%�ַ����ʿռ䣿
                                    delta_x = abs(front_bike(j,3)-back_bike(1,3));
                                    delta_y = abs(front_bike(j,4)-back_bike(1,4));
                                    yy = delta_y.^2/r1 + delta_x.^2/r2; 
                                    if yy<=1 & (front_bike(j,4)-back_bike(1,4)>=0)                                    
                                    overtake_number = overtake_number+1;%���������ۼ�
                                    time_overtake =[time_overtake; t];%����ʱ��
                                    longitudinal_distance =[ longitudinal_distance; abs(a2(1)-a1(1))];%������m
                                    lateral_space =[ lateral_space; max(a4(1)-2,6-a4(1))];%������ʼʱǰ����߽�֮��ļ�����ֵ
                                    lateral_distance =[lateral_distance; abs( front_bike1(j,3)-back_bike1(1,3))];%
                                    if  front_bike(j,3)-back_bike(1,3)>=0
                                        left_msfm = left_msfm+1;
                                    end
                                    tip = 1;
                                    bike_number=[bike_number;i];%��������-��
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
    %�������
    times
    overtake_number
    left_msfm
    xlswrite('��������1',[bike_number,time_overtake,longitudinal_distance,lateral_space,lateral_distance],3);

    %% =======================������������ֲ��Ա�==================
    field_overtake = xlsread('��������1.xls',1);
    sfm_overtake = xlsread('��������1.xls',2);
    msfm_overtake = xlsread('��������1.xls',3);
    % ������Ϣ��������ͼ
    x1 = field_overtake(:,2:5);
    x2 = sfm_overtake(:,2:5);
    x3 = msfm_overtake(:,2:5);
    %x1����ά�ȣ���x2һ��
    xx1 = zeros(length(x2),4);
    xx1(:,:) = NaN;
    xx1(1:length(x1),:) = x1;
    %x3����ά�ȣ���x2һ��
    xx3 = zeros(length(x2),4);
    xx3(:,:) = NaN;
    xx3(1:length(x3),:) = x3;
    X = [xx1 x2 xx3];%�ۺ���������³�����Ϣ����
%  % ������������ͼ
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
% %        saveas(gcf,['��������ʱ������ͼ','.fig']);    
%  LX_FIELD = X(:,4);
%  LX_MYMODEL = X(:,12)+0.8; 
% figure
%     boxplot([LX_FIELD LX_MYMODEL],{'LX_FIELD', 'LX_MYMODEL'})
%     grid on; 
% %        saveas(gcf,['��������ʱ������ͼ','.fig']);

%% ����ʱ������ͼ
    figure
    subplot(1,2,1)
    boxplot([X(:,2) X(:,6) X(:,10)],{'FIELDlongitidual', 'SFMlongitidual','MYMODELlongitidual'})
    grid on;
    subplot(1,2,2)
    boxplot([X(:,3) X(:,7) X(:,11)],{'FIELDclearance', 'SFMclearance','MYMODELclearance'})
    grid on;
    saveas(gcf,['������������ͼ','.fig']);   
    
figure
    boxplot([X(:,4) X(:,8) X(:,12)+0.8],{'FIELDlateral', 'SFMlateral', 'MYMODELlateral'})
    grid on;
   saveas(gcf,['��������ʱ������ͼ','.fig']);

    
%% --------------------ָ��3˲ʱλ�þ���������------------------------
% ����������׼��������ֵ���ƽ���͵�ƽ��ֵ��ƽ������
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
    
      %% ��������������  
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
    %% =====================����Ч��ͼ============================
    data = xlsread('link-ebike','���賬��������Ų��Լ�');
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
%        saveas(gcf,['���������켣',num2str(i),'.fig']);
  saveas(gcf,['���������켣',num2str(i),'.png']);
             hold off;
    end    
    
    %% ======================�켣λ��ͼSx-Sy ======================
    % ��š�ʱ�䡢����λ�á�����λ�á���������ٶȡ���������ٶȡ����������ٶȡ����������ٶ�
    %--------------------ָ��1�켣ͼ------------------------
    %% ��̬����
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
%         saveas(gcf,['�����켣',num2str(i),'.fig']);
 saveas(gcf,['�����켣',num2str(i),'.png']);
    end
    %% ����켣
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
     saveas(gcf, '����켣','fig');
    %% �ٶ�ָ��
    %% --------------------�����ٶȷֲ�-----------------------
    %% ����ƽ���ٶȷֲ����������ʱ��
    speed = zeros(length(Trj_real(:,1)),4);
    [C, ia, ic] = unique(Trj_real(:,1),'rows');%ȡia
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
    [a1,b1]=hist(speed(:,2));%�ٶ�,,,��ͼ���ֶ��ڹ������е����������
    bar(b1,a1/sum(a1));
    [a2,b2]=hist(speed(:,3));%�ٶ�
    bar(b2,a2/sum(a2));
    [a3,b3]=hist(speed(:,4));%�ٶ�
    bar(b3,a3/sum(a3));
    title('The distribution of travel speed');ylabel('Frequency');xlabel('travel speed m/s');
       saveas(gcf, 'The distribution of travel speed','fig');
    %% ----------
    %% ÿ��������� compute
    num = min([length(Trj_real),length(Trj_sfm),length(Trj_MSFMcomfortzone)]);
    t_start = Trj_MSFMcomfortzone(1,2);
    t_end = max(Trj_MSFMcomfortzone(:,2));
    L_REAL = zeros(num,2);%��ࡢƽ���ٶ�
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
                    inter_length_real = sqrt(diff(Trj_real(p_index_real,3)).^2+diff(Trj_real(p_index_real,4)).^2);%�������%real
                    inter_length_sfm = sqrt(diff(Trj_sfm(p_index_sfm,3)).^2+diff(Trj_sfm(p_index_sfm,4)).^2);%�������%SFM
                    inter_length_sfm_fix_passing = sqrt(diff(Trj_MSFMcomfortzone(p_index_sfm_fix_passing,3)).^2+diff(Trj_MSFMcomfortzone(p_index_sfm_fix_passing,4)).^2);%�������%�޸�SFM+FIX+PASSING
                    %ÿ������Ӧ����Χ�����ٶ�
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
  % ���ֱ��ͼ�ֲ�
    L_REAL = L_REAL (1:length(L_REAL),:);
    L_SFM = L_SFM (1:length(L_SFM),:);
    L_SFM_FIX_PASSING = L_SFM_FIX_PASSING (1:length(L_SFM_FIX_PASSING),:);
    L_REAL = L_REAL(L_REAL (:,1)>0.8,:);
    L_SFM = L_SFM(L_SFM (:,1)>0.8,:);
    L_SFM_FIX_PASSING = L_SFM_FIX_PASSING(L_SFM_FIX_PASSING (:,1)>0.8,:);
    %%���ݷ���
    % X = [0.5 2 3.5 5 6.5 8 9.5 11 12.5 14];%���m
    X = [1 2 3 4 5 6 7 8 9 10 ];%���1m
    numgroup = zeros(11,3);%��Ϊ����
    % ��һ��REAL���ڶ���SFM��������SFM_FIX_PASSING
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
%     saveas(gcf,'����ʱ��֡��������','fig');
           
    % ���-�ٶȹ�ϵ�ֲ�   ��ʲô��ȷ��ʵ�������𣿺���û��ʲô���壡����
    L_REAL = L_REAL (1:length(L_REAL),:);
    L_SFM = L_SFM (1:length(L_SFM),:);
    L_SFM_FIX_PASSING = L_SFM_FIX_PASSING (1:length(L_SFM_FIX_PASSING),:);
    %%���ݷ���
    X = [0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5];%���m
    speedgroup = zeros(15,4);%��Ϊ����
    % ��һ��REAL���ڶ���SFM��������SFM_FIX��������SFM_FIX_PASSING
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
    title '�����ٶȾ�ֵ����Ĺ�ϵ';
    set(gca,'XTicklabel',{'0.5' '1.5' '2.5' '3.5' '4.5' '5.5' '6.5' '7.5' '8.5' '9.5' '10.5' '11.5' '12.5' '13.5' '14.5'});
    h = legend('FIELD','SFM','MSFM');
%     h = legend('FIELD','MYMODEL');    
%     set(h, 'Box', 'off');
%     set(h,'Interpreter','none')
    xlabel('�������/m'); ylabel('ƽ���ٶ�m/s');
    % % saveas(gcf,'����ʱ��֡��������','jpg');
    axis([0.5 13 4 10]);
  grid;
       saveas(gcf,'�����ٶȾ�ֵ����Ĺ�ϵ','fig');  

% % ���������ֲ�   
% figure
% xbins = 0:5:30;
% % histogram(LY_FIELD,xbins);
% h = histogram(LY_FIELD,'Normalization','probability');
% h.BinEdges = [0:5:30];
% hold on
% % histogram(LY_MYMODEL,xbins);
% h = histogram(LY_MYMODEL,'Normalization','probability');
% h.BinEdges = [0:5:30];
% %  saveas(gcf,['������ʼ�׶�������ֲ�ͼ','.fig']);
% figure
% xbins = 1:0.5:5;
% % histogram(CLEARENCE_FIELD,xbins);
% h = histogram(CLEARENCE_FIELD,'Normalization','probability');
% h.BinEdges = [1:0.5:5];
% hold on
% % histogram(CLEATENCE_MYMODEL,xbins);
% h = histogram(CLEATENCE_MYMODEL,'Normalization','probability');
% h.BinEdges = [1:0.5:5];
% %   saveas(gcf,['������ʼ�׶κ��򾻿շֲ�ͼ','.fig']); 
% figure
% xbins = 0:0.1:3;
% % histogram(LX_FIELD,xbins);
% h = histogram(LX_FIELD,'Normalization','probability');
% h.BinEdges = xbins;
% hold on
% % histogram(LX_MYMODEL,xbins);
% h = histogram(LX_MYMODEL,'Normalization','probability');
% h.BinEdges = xbins;
% %    saveas(gcf,['�������к�����ʱ�̷ֲ�ͼ','.fig']);
%    