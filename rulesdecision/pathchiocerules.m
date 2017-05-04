%% pathchoicesimulation 
function [tagA,tagB,tag_rang,tag_emergency,pathchoice_dataset]= pathchiocerules(subject_vehicle,vehicle,LX,LY,LYB,leadvehicle,TLB)
tagA =0; tagB = 0; tag_rang=0; tag_emergency = 0;%���� 
vehicle1 = vehicle(vehicle(:,1)~=subject_vehicle(1,1),:);
pathchoice_dataset = []; pathchoice_dataset_final = [];
    lead_vehicle_number = leadvehicle(1,1);%  
    longitudinal = leadvehicle(1,4)-subject_vehicle(1,4);
    speed_subject =  subject_vehicle(1,6);
    speed_lead = leadvehicle(1,6);
    delta_speed_sp = -(speed_subject-speed_lead);
    distance_sp = sqrt(sum(((subject_vehicle(1,3:4)-leadvehicle(1,3:4)).^2).'));   
    lateral_sp = subject_vehicle(1,3)-leadvehicle(1,3);
         if  lateral_sp>0
              lateral_spl = lateral_sp;
              lateral_spr = -lateral_spl;              
          else
              lateral_spr = -lateral_sp;
              lateral_spl = -lateral_spr;                 
         end         
           clearance_pl = abs(6-leadvehicle(1,3));
           clearance_pr = abs(2-leadvehicle(1,3));   
           LMAX = abs(subject_vehicle(1,4)-leadvehicle(1,4));
          leadaside_vehicles = vehicle1(vehicle1(:,4)-leadvehicle(1,4)>-LMAX&vehicle1(:,4)-leadvehicle(1,4)<=2,:);%ǰ���෽���г���
          if ~isempty(leadaside_vehicles)
          leadleft_vehicle = leadaside_vehicles(leadaside_vehicles(:,3)-leadvehicle(1,3)>0,:);%ǰ�����ǰ������
          if ~isempty(leadleft_vehicle)
          leadleft_vehicle = leadleft_vehicle(leadleft_vehicle(:,3)==min(leadleft_vehicle(:,3)),:);%·����ǰ�����ǰ������
          leadleft_vehicle = leadleft_vehicle(leadleft_vehicle(:,4)==min(leadleft_vehicle(:,4)),:);%·����ǰ�����ǰ������  =============        
          %������ǰ����೵���ٶȲ�
             delta_speed_spl = -(speed_subject-leadleft_vehicle(1,6)); 
          %������ǰ����೵������
             distance_spl= sqrt(sum(((subject_vehicle(1,3:4)-leadleft_vehicle(1,3:4)).^2).'));
          %ǰ��������೵���ĺ�����
          lateral_pl = abs(leadleft_vehicle(1,3)-leadvehicle(1,3));             
          else
          %������ǰ����೵���ٶȲ� 
              delta_speed_spl = speed_subject-0;  
          %������ǰ����೵������
             distance_spl= 25;%����೵����ȡһ����ֵ
          %ǰ��������೵���ĺ�����
          lateral_pl = min(abs(6-leadvehicle(:,3)));              
          end
          %Ѱ��ǰ���Ҳ෽����          
          leadright_vehicle = leadaside_vehicles(leadaside_vehicles(:,3)-leadvehicle(1,3)<0,:);%ǰ���Ҳ�ǰ������
          if ~isempty(leadright_vehicle)
          leadright_vehicle = leadright_vehicle(leadright_vehicle(:,3)==min(leadright_vehicle(:,3)),:);%·����ǰ���Ҳ�ǰ������
          leadright_vehicle = leadright_vehicle(leadright_vehicle(:,4)==min(leadright_vehicle(:,4)),:);%·����ǰ���Ҳ�ǰ������          
          %������ǰ���Ҳ೵���ٶȲ�
             delta_speed_spr = -(speed_subject-leadright_vehicle(1,6)); 
          %������ǰ���Ҳ೵������
             distance_spr= sqrt(sum(((subject_vehicle(1,3:4)-leadright_vehicle(1,3:4)).^2).'));
          %ǰ�������Ҳ೵���ĺ�����
          lateral_pr = abs(leadright_vehicle(1,3)-leadvehicle(1,3));     
          else
          %������ǰ���Ҳ೵���ٶȲ� 
              delta_speed_spr = speed_subject-0;  
          %������ǰ���Ҳ೵������
             distance_spr= 25;%���Ҳ೵����ȡһ����ֵ
          %ǰ�������Ҳ�߽�ĺ�����
          lateral_pr = min(abs(2-leadvehicle(:,3)));   
          end   
          else
          %������ǰ����೵���ٶȲ� 
              delta_speed_spl = speed_subject-0;  
          %������ǰ����೵������
             distance_spl= 25;%����೵����ȡһ����ֵ
          %ǰ��������೵���ĺ�����
          lateral_pl = min(abs(6-leadvehicle(:,3))); 
          %������ǰ���Ҳ೵���ٶȲ� 
              delta_speed_spr = speed_subject-0;  
          %������ǰ���Ҳ೵������
             distance_spr= 25;%���Ҳ೵����ȡһ����ֵ
          %ǰ�������Ҳ೵���ĺ�����
          lateral_pr = min(abs(2-leadvehicle(:,3)));           
          end   
  %Ѱ�ұ����෽����(����������г���) 
           subjectaside_vehicles = vehicle1(vehicle1(:,4)-subject_vehicle(1,4)>=-2&vehicle1(:,4)-subject_vehicle(1,4)<=0,:);%�����෽���г���
           if ~isempty(subjectaside_vehicles)
              subjectleft_vehicle = subjectaside_vehicles(subjectaside_vehicles(:,3)-subject_vehicle(1,3)>0,:);
              if ~isnan(subjectleft_vehicle)%������೵��
              subjectleft_vehicle = subjectleft_vehicle(subjectleft_vehicle(:,3)==min(subjectleft_vehicle(:,3)),:);
              subjectleft_vehicle = subjectleft_vehicle(subjectleft_vehicle(:,4)==max(subjectleft_vehicle(:,4)),:);              
              %����������೵�������󷽣��ٶȲ�
              delta_speed_sl = speed_subject-subjectleft_vehicle(1,6);
              %����������೵������
             distance_sl= sqrt(sum(((subject_vehicle(1,3:4)-subjectleft_vehicle(1,3:4)).^2).'));
              else
              %����������೵�������󷽣��ٶȲ�
              delta_speed_sl = speed_subject-0;
              %����������೵������
              distance_sl= 10;%����������೵������ȡ����ֵ
              end
             %Ѱ���Ҳ೵��
              subjectright_vehicles = subjectaside_vehicles(subjectaside_vehicles(:,3)-subject_vehicle(1,3)<0,:);
              if ~isnan(subjectright_vehicles)%�����Ҳ೵��
              subjectright_vehicle = subjectright_vehicles(subjectright_vehicles(:,3)==min(subjectright_vehicles(:,3)),:);
              subjectright_vehicle = subjectright_vehicles(subjectright_vehicles(:,4)==max(subjectright_vehicles(:,4)),:);              
              %���������Ҳ೵�����Ҳ�󷽣��ٶȲ�
              delta_speed_sr = speed_subject-subjectright_vehicle(1,6);
              %���������Ҳ೵������
             distance_sr= sqrt(sum(((subject_vehicle(1,3:4)-subjectright_vehicle(1,3:4)).^2).'));         
              else
              %���������Ҳ೵�����Ҳ�󷽣��ٶȲ�
              delta_speed_sr = speed_subject-0;
              %���������Ҳ೵������
              distance_sr= 10;%���������Ҳ೵������ȡ����ֵ
              end   
           else
              %����������೵�������󷽣��ٶȲ�
              delta_speed_sl = speed_subject-0;
              %����������೵������
              distance_sl= 10;%����������೵������ȡ����ֵ
              %���������Ҳ೵�����Ҳ�󷽣��ٶȲ�
              delta_speed_sr = speed_subject-0;
              %���������Ҳ೵������
              distance_sr= 10;%���������Ҳ೵������ȡ����ֵ              
           end   
           tag=0;
    pathchoice_dataset = [subject_vehicle(1,1), lead_vehicle_number,subject_vehicle(1,2),tag...
                         delta_speed_sp,distance_sp, lateral_spl ,lateral_spr,clearance_pl,clearance_pr,delta_speed_spl,...
                           distance_spl,lateral_pl,delta_speed_spr,distance_spr,lateral_pr,delta_speed_sl,distance_sl,delta_speed_sr,distance_sr];%·��ѡ����Ϣ
   pathchoice_dataset_final = [pathchoice_dataset_final;pathchoice_dataset];        
   
lateral_distance = max([clearance_pl;clearance_pr;lateral_pr;lateral_pl]);
   if  (longitudinal>LYB&longitudinal<LY)%A�㴦 �ٶȺ�λ��΢�����ӽ��׶�       
%��������   
      if (delta_speed_sp<=0 & lateral_distance>=LX)|(delta_speed_sp<=-2 & lateral_distance>=max(0.8*LX,1))%overtaking conditions forwards
         if (delta_speed_sl>=0&distance_sl>1.2)|(delta_speed_sr>=0&distance_sr>1.2)%overtaking conditions surrounding
%overtaking direction
            if  (max([clearance_pl;lateral_pl])>max([clearance_pr;lateral_pr]))&(delta_speed_sl>=0&distance_sl>1.2)% turning left
                tagA= 1;
            else if (max([clearance_pl;lateral_pl])<=max([clearance_pr;lateral_pr]))&(delta_speed_sr>=0&distance_sr>1.2)
                   tagA = 3; 
                else 
                  tagA = 2;
                end
            end 
         else
             tagA =2;
         end   
      else
          tagA =2;
      end
   end  
  
   if (longitudinal<=LYB&longitudinal>=3) %B�� ����ִ��      
      if (delta_speed_sp<=0 & lateral_distance>=LX)|(delta_speed_sp<=-2 & lateral_distance>=max(0.8*LX,1))%overtaking conditions forwards    
         if (delta_speed_sl>=0&distance_sl>1.2)|(delta_speed_sr>=0&distance_sr>1.2)%overtaking conditions surrounding  
            if  TLB(1,1)==1
                tagB = 1;    
            else if TLB(1,2)==1
                  tagB = 3;
                else if (max([clearance_pl;lateral_pl])>max([clearance_pr;lateral_pr]))&(delta_speed_sl>=0&distance_sl>1.2)% turning left
                        tagB= 1;
                    else if (max([clearance_pl;lateral_pl])<=max([clearance_pr;lateral_pr]))&(delta_speed_sr>=0&distance_sr>1.2)
                         tagB= 3;
                        else
                         tagB = 2;     
                        end
                    end   
                end
            end 
           tagB = 2;  
         end
      else
          tagB = 2;     
      end    
   end  
     
    if longitudinal<=max([3;LX+1.8]) %ǰ���������������%4
        tag_rang = 1;
    end
     if longitudinal<1.5 & abs(leadvehicle(1,3)-subject_vehicle(1,3))<0.5  %C�� �������
       tag_emergency = 1;
     end
end