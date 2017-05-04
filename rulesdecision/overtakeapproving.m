%% �˶�ִ��ģ�� 
function       [newacceleration] = overtakeapproving(subject_vehicle,leadvehicle)
    newspeed = zeros(1,2);
    newacceleration = zeros(1,2);
    %the distance between subject_vehicle and leadvehicle      
   distance = subject_vehicle(1,4)-leadvehicle(1,4);
   
   if distance<=-6 %�����ӽ��׶�
    lateral_distance = abs(subject_vehicle(1,3)-leadvehicle(1,3));   
     if  lateral_distance<0.6
      delta_speed_lateral = -0.0032*distance-0.21;
      delta_speed_longit = -0.21*distance+0.72;
      newspeed(1,1) = leadvehicle(1,5)+delta_speed_lateral;
      newspeed(1,2) = leadvehicle(1,6)+delta_speed_longit;
      newacceleration = (newspeed-subject_vehicle(1,5:6))./0.12; 
     end     
   else   %����ƫ�ƺͺ��򱣳ֽ׶�
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
           if distance<11  %��ʵ�ǲ������õģ���û�з���ԭ�г����Ļ��ƣ���Ҫ�������¼�����ĳ�������ڲ�ͬʱ�̣��ǲ���һ��
   delta_speed_lateral = -0.012*distance+0.23;
   delta_speed_longit = 0.11*distance+1.5;
   newspeed(1,1) = leadvehicle(1,5)+delta_speed_lateral;
   newspeed(1,2) = leadvehicle(1,6)+delta_speed_longit;
   newacceleration = (newspeed-subject_vehicle(1,5:6))./0.12;               
           end
       end  
   end
end
