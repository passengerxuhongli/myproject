%期望加速度的形式表征期望速度运行的阶段
function [a_driving] = drivingforce(vehicle)
    tao = 3.87;
 % long-term planning
    %期望目的地方向向量    
    r_d = sqrt( (vehicle(:,9) - vehicle(:,3)).^2 + (vehicle(:,10) - vehicle(:,4)).^2 );
    exp_d = [(vehicle(:,9) - vehicle(:,3)) ./ r_d, (vehicle(:,10) - vehicle(:,4)) ./ r_d];
   %驱动加速度   
    a_driving(:,1) = ( vehicle(:,11).*exp_d(:,1) - vehicle(:,5) ) ./ tao;
    a_driving(:,2) = ( vehicle(:,11).*exp_d(:,2) - vehicle(:,6) ) ./ tao;
end