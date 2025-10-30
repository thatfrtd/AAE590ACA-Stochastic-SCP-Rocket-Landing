
%% Covariance
nx = 4;
nu = 2;
B = [zeros(2); eye(2)];
c = zeros(nx,1);

A = @(t, x) A_func(t, x, mu_sun);

i = 1;
Pk_his = zeros(nx,nx,length(xr));
Pk_his(:,:,1) = P_0;
while i <= (length(xr) - 1)
    
    %[dc, u, A] = lt_ot_con_2d(tr(i),xr(:,i),mu_sun);
    
    %[dc, A] = lt_ot_con_2d_brown(t,c,mu,dt,G, w_his(ii*2-1:ii*2,i) , t_his(i), u_his(:,i));
 
    %u = @(t) interp1([tr(i), tr(i+1)], [u_his(:, i), u_his(:, i + 1)], t);
    
    y02 = eye(nx);
    y03 = zeros(nx,nu);
    y04 = zeros(nx,1);
    y05 = zeros(nx,nx);
    Y0 = [xr(1:nx,i); y02(:); y03(:); y04(:); y05(:)];
    
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [t, traj] = ode45(@(t, Y) EoMwithDiscreteMatrixG(t, Y, u_his(:, i), A, B, c, nx, G), [tr(i), tr(i+1)], Y0, options);
    
    Yf = traj(end,:);
    xf = Yf(1:nx);
    
    Ak = reshape(Yf(nx+1:nx+nx^2),nx,nx);
    Bk = Ak*reshape(Yf(nx+nx^2+1:nx+nx^2+nx*nu),nx,nu);
    ck = Ak*(Yf(nx+nx^2+nx*nu+1:nx+nx^2+nx*nu+nx))';

    sumk = Ak*reshape(Yf(:, (nx*(nx+1) + nx*nu + nx)+1:(nx*(nx+1) +nx*nu + nx)+nx*nx), nx, nx)*Ak';
    Gk = chol(sumk, 'lower');
    Pk = Ak*Pk_his(:, :, i)*Ak'+Gk*Gk';
    Pk_his(:, :, i+1) = Pk;

    i = i + 1;
end

dR1 = R1 - xr(1,:);
dR2 = R2 - xr(2,:);
dV1 = V1 - xr(3,:);
dV2 = V2 - xr(4,:);

diagonalP = [];
for i=1:length(Pk_his)
    diagonalP(:,i) = diag(Pk_his(i));
end
three_sig_p = 3*sqrt(diagonalP);
function [dc, u, A] = lt_ot_con_2d(t,c,mu)
    x1 = c(1);
    x2 = c(2);
    x3 = c(3);
    x4 = c(4);
    r = [x1 x2]';
    v = [x3 x4]';
    a = norm(r);
    l = [c(5) c(6) c(7) c(8)]';
    B = [zeros(2); eye(2)];
    
    %u = (-1/2)*B'*l; % for min energy

    %{
    %for min time
    p = -B'*l; % for min time
    umax = 0.1;
    u = (p/norm(p))*umax;
    %}

    %{
    uhat = p/norm(p);
    if norm(p) >1
        Gam = umax;
    elseif norm(p) <1
        Gam = 0;
    end
    u = Gam*uhat;
    %}

    %for min fuel
    
    p = -B'*l;
    uhat = p/norm(p);
    umax = 0.1;
    if 1-norm(p)<0
        u = umax*uhat;
    elseif 1-norm(p)>0
        u = 0*uhat;
    end
    
    dc = [v; -mu*r/(a^3)]+B*u;
    A = [zeros(2) eye(2); (eye(2) - 3 * (r * r') / a^2) * (-mu/(a^3)) zeros(2)];
    %dc = [v; -r/(a^(3/2))]+B*u;
    

    l_dot = -[zeros(2) eye(2); -(eye(2)/(a^3) - (3/(a^5))*(r*r')) zeros(2)]'*l;

    dc = [dc; l_dot];
    %dc = [dc; reshape(KDot, [16 1])];
return
end

function [A] = A_func(t, x, mu)
    r = x(1:2);
    a = norm(r);
    A = [zeros(2) eye(2); (eye(2) - 3 * (r * r') / a^2) * (-mu/(a^3)) zeros(2)];
end

function [dx, phi] = EoMwithDiscreteMatrixG(t, Y, uk,A,B,c, n, G)
    x = Y(1:n);
    phi = reshape(Y(n+1:n+n^2),n,n);
    dc2 = A(t,x)*phi;
    dc3 = phi\B;
    dc4 = phi\c;
    dc5 = (phi\G)*(phi\G)';
    dx = [A(t,x)*x+B*uk+c; dc2(:); dc3(:); dc4(:); dc5(:)];

end