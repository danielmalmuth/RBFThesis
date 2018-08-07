function q_new = rk4_adv(delta_t,u,D,q_init,target_time)

q_new = q_init;

for i = 1:(target_time/delta_t)
    k1 = (delta_t*u).*D*q_new;
    k2 = (delta_t*u).*D*(q_new + k1./2);
    k3 = (delta_t*u).*D*(q_new + k2./2);
    k4 = (delta_t*u).*D*(q_new + k3);
    q_new = q_new + (k1 + 2.*k2 + 2.*k3 + k4)./6;
end

end