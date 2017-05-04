%% 计算舒适空间的长短半轴，半径等信息
function  [comfort_lateral,interdistance,LX,LY,LYB,leadvehicle] = comfortzonefeature(subject_vehicle,vehicle)
        Q = sum(vehicle(:,4)-subject_vehicle(1,4)>=-5 & vehicle(:,4)-subject_vehicle(1,4)<=5);
        density = Q/(4*10);
        interdistance = min(0.156*density^(-0.18746),1.6);
        comfort_lateral = min(0.8+interdistance,1.6);
        lx = 0.8-0.5*subject_vehicle(1,11)*sin(pi-2*acos(9.8*1.0*tan(pi*15/180)/(subject_vehicle(1,11).^2)));
        LX = lx+comfort_lateral;
%前车
       vehicle1 =  vehicle(vehicle(:,1)~=subject_vehicle(1,1),:);
       leadvehicle = vehicle1(abs(vehicle1(:,3)-subject_vehicle(1,3))<=LX...
           & vehicle1(:,4)-subject_vehicle(1,4)>=0,:);
       v_subjectvehicle = subject_vehicle(1,6);
       if ~isempty(leadvehicle)
           leadvehicle = leadvehicle(leadvehicle(:,4)==min(leadvehicle(:,4)),:);
           v_leadvehicle = leadvehicle(1,6);
           Q1 = sum(abs(vehicle(:,4)-leadvehicle(1,4))<=2.5);
           k_leadvehicle = Q1/(4*5);
           Q2 = sum(abs(vehicle(:,4)-subject_vehicle(1,4))<=2.5);
           k_subjectvehicle = Q2/(4*5);
       end
       if  ~isempty(leadvehicle) 
            v_leadvehicle = leadvehicle(1,6);
           if v_subjectvehicle>=v_leadvehicle
              LY = abs(15.3911+0.5*v_subjectvehicle+0.1784*v_subjectvehicle.^2-...
                        0.4234*v_leadvehicle.^2-56.1606*k_leadvehicle);
              LYB = abs(8.628+0.613*v_subjectvehicle-0.42*v_leadvehicle-31.628*k_subjectvehicle);      
           else 
              LY =  v_subjectvehicle*0.5+v_subjectvehicle^2/(2*2);  
              LYB = max(0.5*LY,3);
           end
      else 
            LY =  v_subjectvehicle*0.5+v_subjectvehicle^2/(2*2);   
            LYB = 0.5*LY;
       end
       %% 侵犯舒适空间的前车
       delta_x = abs(vehicle1(:,3)-subject_vehicle(1,3));
       delta_y = abs(vehicle1(:,4)-subject_vehicle(1,4));
       yy = delta_y.^2/LY + delta_x.^2/LX;
       leadvehicle = vehicle1(yy<=1 & vehicle1(:,4)-subject_vehicle(1,4)>=0,:);
       
       
end