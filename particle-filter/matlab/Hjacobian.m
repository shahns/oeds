    function H = Hjacobian(x,x0)
% H = Hjacobian(x,hparm) gives the 3x8 Jacobian matrix of the measurement
% function evaluated at the state x. x0 is the position of the beacon.

H = zeros(3,8);
H(1,1) =  (180/pi)*(x0(2)-x(2))./((x0(1)-x(1)).^2 + (x0(2)-x(2)).^2);
H(1,2) = -(180/pi)*(x0(1)-x(1))./((x0(1)-x(1)).^2 + (x0(2)-x(2)).^2);

s = sqrt(x(3)^2+x(4)^2);
if s==0
    
    s= 10*eps;
    x(3) = 10*eps;
    x(4) = 10*eps;
end
s = sqrt(x(3)^2+x(4)^2);

H(2,3) = x(3)/s;
H(2,4) = x(4)/s;
H(3,3) = -(180/pi)*x(4)/s^2;
H(3,4) =  (180/pi)*x(3)/s^2;
