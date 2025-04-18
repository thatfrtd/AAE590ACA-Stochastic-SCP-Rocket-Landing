function xdot = SymDynamics3DoF(t,in2,in3,mass,L,I)
%SymDynamics3DoF
%    XDOT = SymDynamics3DoF(T,IN2,IN3,MASS,L,I)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    13-Apr-2025 02:57:54

theta1 = in2(5,:);
thrust1 = in3(1,:);
thrust2 = in3(2,:);
v1 = in2(3,:);
v2 = in2(4,:);
w1 = in2(6,:);
t2 = cos(theta1);
t3 = sin(theta1);
t4 = 1.0./mass;
xdot = [v1;v2;t4.*(t2.*thrust1-t3.*thrust2);t4.*(t2.*thrust2+t3.*thrust1)-3.7114e-3;w1;-(L.*thrust2)./I];
end
