%% ��������ʱ��������ٶȣ�λ��
function [start_state] = initdata(bike)
   index = find(~isnan(bike(:,1)));
   start_state = zeros(length(index),11);%���ɳ����ĳ�ʼ״̬
     for i = 1:length(index)
         start_state(i,1) = bike(index(i),1);   
         start_state(i,2:8) = bike(index(i)+1,2:8);
         start_state(i,9:10) = [normrnd(3.22,1.72),120];%normrnd(3.22,0.72)
         while start_state(i,9)<2|start_state(i,9)>6
           start_state(i,9:10) = [normrnd(3.22,1.72),120];%normrnd(3.22,0.72)
         end
         start_state(i,11) = normrnd(7.1,1.93);
     end
     t0 = min(start_state(:,2));%�����ʼʱ��
     start_state(:,2) = (start_state(:,2)-t0)*24*3600;
end