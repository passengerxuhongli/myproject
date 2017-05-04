%% 运动执行模块 
function       [newacceleration] = overtakeapproving(subject_vehicle,leadvehicle)
    newspeed = zeros(1,2);
    newacceleration = zeros(1,2);
    %the distance between subject_vehicle and leadvehicle      
   distance = subject_vehicle(1,4)-leadvehicle(1,4);
   
   if distance<=-6 %超车接近阶段
    lateral_distance = abs(subject_vehicle(1,3)-leadvehicle(1,3));   
     if  lateral_distance<0.6
      delta_speed_lateral = -0.0032*distance-0.21;
      delta_speed_longit = -0.21*distance+0.72;
      newspeed(1,1) = leadvehicle(1,5)+delta_speed_lateral;
      newspeed(1,2) = leadvehicle(1,6)+delta_speed_longit;
      newacceleration = (newspeed-subject_vehicle(1,5:6))./0.12; 
     end     
   else   %超车偏移和横向保持阶段
       if  distance<=3%0
       lateral_distance = abs(subject_vehicle(1,3)-leadvehicle(1,3));   
        if lateral_distance<1.5
        delta_speed_lateral = 0.038*distance-0.002;
        delta_speed_longit = 0.11*distance+1.5;
        newspeed(1,1) = leadvehicle(1,5)+delta_speed_lateral;
        newspeed(1,2) = leadvehicle(1,6)+delta_speed_longit;
        newacceleration = (newspeed-subject_vehicle(1,5:6))./0.12;         
        end   
       else 
           if distance<11  %其实是不起作用的，即没有返回原有车道的机制，若要考虑需记录超车的车辆编号在不同时刻，是不是一致
   delta_speed_lateral = -0.012*distance+0.23;
   delta_speed_longit = 0.11*distance+1.5;
   newspeed(1,1) = leadvehicle(1,5)+delta_speed_lateral;
   newspeed(1,2) = leadvehicle(1,6)+delta_speed_longit;
   newacceleration = (newspeed-subject_vehicle(1,5:6))./0.12;               
           end
       end  
   end
end
