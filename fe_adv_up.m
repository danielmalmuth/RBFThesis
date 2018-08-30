function q_new = fe_adv_up(delta_t,u_func,q_init,target_time,h)

q_new = q_init;

for i = 1:(target_time/delta_t)
    u = u_func((i-1)*delta_t);
    Dq = up_adv_diff(u*delta_t,q_new,h);
    q_new = q_new + (delta_t*u).*Dq;
    plot(q_new,'k-')
    ylim([0 3])
    pause(0.1);
end

end