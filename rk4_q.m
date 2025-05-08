function position = rk4(x, u, c, t)

    function q_integrated = q_deriv(q, dq, t)
        q_integrated = q + 0.5 * t * [dq(1) * q(1) - dq(2) * q(2) - dq(3) * q(3) - dq(4) * q(4)
                                           dq(1) * q(2) + dq(2) * q(1) - dq(3) * q(4) + dq(4) * q(3)
                                           dq(1) * q(3) + dq(2) * q(4) + dq(3) * q(1) - dq(4) * q(2)
                                           dq(1) * q(4) - dq(2) * q(3) + dq(3) * q(2) + dq(4) * q(1)].';
        q_integrated = q_integrated / norm(q_integrated);
    end    

    timestep = t / 2;
    

    k0 = SymDynamics6DoF_q(x, u, c)
    est1 = [ timestep .* k0(1:6).' + x(1:6) q_deriv(x(7:10), k0(7:10).', timestep) timestep .* k0(11:13).' + x(11:13) ];

    k1 = SymDynamics6DoF_q(est1, u, c)
    est2 = [ timestep * k1(1:6).' + x(1:6) q_deriv(x(7:10), k1(7:10).', timestep) timestep * k1(11:13).' + x(11:13) ];
    
    k2 = SymDynamics6DoF_q(est2, u, c)
    position = [ t * k2(1:6).' + x(1:6) q_deriv(x(7:10), k2(7:10).', t) t * k2(11:13).' + x(11:13) ];

    k3 = SymDynamics6DoF_q(position, u, c)

    position = x + (t/6) * (k0.' + 2 * k1.' + 2 * k2.' + k3.');
    
    

end