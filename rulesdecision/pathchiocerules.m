%% pathchoicesimulation 
function [tagA,tagB,tag_rang,tag_emergency,pathchoice_dataset]= pathchiocerules(subject_vehicle,vehicle,LX,LY,LYB,leadvehicle,TLB)
tagA =0; tagB = 0; tag_rang=0; tag_emergency = 0;%归零 
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
          leadaside_vehicles = vehicle1(vehicle1(:,4)-leadvehicle(1,4)>-LMAX&vehicle1(:,4)-leadvehicle(1,4)<=2,:);%前车侧方所有车辆
          if ~isempty(leadaside_vehicles)
          leadleft_vehicle = leadaside_vehicles(leadaside_vehicles(:,3)-leadvehicle(1,3)>0,:);%前车左侧前方车辆
          if ~isempty(leadleft_vehicle)
          leadleft_vehicle = leadleft_vehicle(leadleft_vehicle(:,3)==min(leadleft_vehicle(:,3)),:);%路径上前车左侧前方车辆
          leadleft_vehicle = leadleft_vehicle(leadleft_vehicle(:,4)==min(leadleft_vehicle(:,4)),:);%路径上前车左侧前方车辆  =============        
          %本车与前车左侧车辆速度差
             delta_speed_spl = -(speed_subject-leadleft_vehicle(1,6)); 
          %本车与前车左侧车辆距离
             distance_spl= sqrt(sum(((subject_vehicle(1,3:4)-leadleft_vehicle(1,3:4)).^2).'));
          %前车与其左侧车辆的横向间距
          lateral_pl = abs(leadleft_vehicle(1,3)-leadvehicle(1,3));             
          else
          %本车与前车左侧车辆速度差 
              delta_speed_spl = speed_subject-0;  
          %本车与前车左侧车辆距离
             distance_spl= 25;%无左侧车辆则取一个极值
          %前车与其左侧车辆的横向间距
          lateral_pl = min(abs(6-leadvehicle(:,3)));              
          end
          %寻找前车右侧方车辆          
          leadright_vehicle = leadaside_vehicles(leadaside_vehicles(:,3)-leadvehicle(1,3)<0,:);%前车右侧前方车辆
          if ~isempty(leadright_vehicle)
          leadright_vehicle = leadright_vehicle(leadright_vehicle(:,3)==min(leadright_vehicle(:,3)),:);%路径上前车右侧前方车辆
          leadright_vehicle = leadright_vehicle(leadright_vehicle(:,4)==min(leadright_vehicle(:,4)),:);%路径上前车右侧前方车辆          
          %本车与前车右侧车辆速度差
             delta_speed_spr = -(speed_subject-leadright_vehicle(1,6)); 
          %本车与前车右侧车辆距离
             distance_spr= sqrt(sum(((subject_vehicle(1,3:4)-leadright_vehicle(1,3:4)).^2).'));
          %前车与其右侧车辆的横向间距
          lateral_pr = abs(leadright_vehicle(1,3)-leadvehicle(1,3));     
          else
          %本车与前车右侧车辆速度差 
              delta_speed_spr = speed_subject-0;  
          %本车与前车右侧车辆距离
             distance_spr= 25;%无右侧车辆则取一个极值
          %前车与其右侧边界的横向间距
          lateral_pr = min(abs(2-leadvehicle(:,3)));   
          end   
          else
          %本车与前车左侧车辆速度差 
              delta_speed_spl = speed_subject-0;  
          %本车与前车左侧车辆距离
             distance_spl= 25;%无左侧车辆则取一个极值
          %前车与其左侧车辆的横向间距
          lateral_pl = min(abs(6-leadvehicle(:,3))); 
          %本车与前车右侧车辆速度差 
              delta_speed_spr = speed_subject-0;  
          %本车与前车右侧车辆距离
             distance_spr= 25;%无右侧车辆则取一个极值
          %前车与其右侧车辆的横向间距
          lateral_pr = min(abs(2-leadvehicle(:,3)));           
          end   
  %寻找本车侧方车辆(本车侧后方所有车辆) 
           subjectaside_vehicles = vehicle1(vehicle1(:,4)-subject_vehicle(1,4)>=-2&vehicle1(:,4)-subject_vehicle(1,4)<=0,:);%本车侧方所有车辆
           if ~isempty(subjectaside_vehicles)
              subjectleft_vehicle = subjectaside_vehicles(subjectaside_vehicles(:,3)-subject_vehicle(1,3)>0,:);
              if ~isnan(subjectleft_vehicle)%存在左侧车辆
              subjectleft_vehicle = subjectleft_vehicle(subjectleft_vehicle(:,3)==min(subjectleft_vehicle(:,3)),:);
              subjectleft_vehicle = subjectleft_vehicle(subjectleft_vehicle(:,4)==max(subjectleft_vehicle(:,4)),:);              
              %本车与其左侧车辆（左侧后方）速度差
              delta_speed_sl = speed_subject-subjectleft_vehicle(1,6);
              %本车与其左侧车辆距离
             distance_sl= sqrt(sum(((subject_vehicle(1,3:4)-subjectleft_vehicle(1,3:4)).^2).'));
              else
              %本车与其左侧车辆（左侧后方）速度差
              delta_speed_sl = speed_subject-0;
              %本车与其左侧车辆距离
              distance_sl= 10;%若不存在左侧车辆，则取个极值
              end
             %寻找右侧车辆
              subjectright_vehicles = subjectaside_vehicles(subjectaside_vehicles(:,3)-subject_vehicle(1,3)<0,:);
              if ~isnan(subjectright_vehicles)%存在右侧车辆
              subjectright_vehicle = subjectright_vehicles(subjectright_vehicles(:,3)==min(subjectright_vehicles(:,3)),:);
              subjectright_vehicle = subjectright_vehicles(subjectright_vehicles(:,4)==max(subjectright_vehicles(:,4)),:);              
              %本车与其右侧车辆（右侧后方）速度差
              delta_speed_sr = speed_subject-subjectright_vehicle(1,6);
              %本车与其右侧车辆距离
             distance_sr= sqrt(sum(((subject_vehicle(1,3:4)-subjectright_vehicle(1,3:4)).^2).'));         
              else
              %本车与其右侧车辆（右侧后方）速度差
              delta_speed_sr = speed_subject-0;
              %本车与其右侧车辆距离
              distance_sr= 10;%若不存在右侧车辆，则取个极值
              end   
           else
              %本车与其左侧车辆（左侧后方）速度差
              delta_speed_sl = speed_subject-0;
              %本车与其左侧车辆距离
              distance_sl= 10;%若不存在左侧车辆，则取个极值
              %本车与其右侧车辆（右侧后方）速度差
              delta_speed_sr = speed_subject-0;
              %本车与其右侧车辆距离
              distance_sr= 10;%若不存在右侧车辆，则取个极值              
           end   
           tag=0;
    pathchoice_dataset = [subject_vehicle(1,1), lead_vehicle_number,subject_vehicle(1,2),tag...
                         delta_speed_sp,distance_sp, lateral_spl ,lateral_spr,clearance_pl,clearance_pr,delta_speed_spl,...
                           distance_spl,lateral_pl,delta_speed_spr,distance_spr,lateral_pr,delta_speed_sl,distance_sl,delta_speed_sr,distance_sr];%路径选择信息
   pathchoice_dataset_final = [pathchoice_dataset_final;pathchoice_dataset];        
   
lateral_distance = max([clearance_pl;clearance_pr;lateral_pr;lateral_pl]);
   if  (longitudinal>LYB&longitudinal<LY)%A点处 速度和位置微调，接近阶段       
%超车规则   
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
  
   if (longitudinal<=LYB&longitudinal>=3) %B点 超车执行      
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
     
    if longitudinal<=max([3;LX+1.8]) %前车横向调整，让行%4
        tag_rang = 1;
    end
     if longitudinal<1.5 & abs(leadvehicle(1,3)-subject_vehicle(1,3))<0.5  %C点 紧急情况
       tag_emergency = 1;
     end
end