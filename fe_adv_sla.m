function q_new = fe_adv_sla(delta_t,u_func,x,q_init,target_time)

q_new = q_init;

for i = 1:(target_time/delta_t)
    u = u_func((i-1)*delta_t);
    x_new = x + -u*delta_t; % new set of x values
    q_new = interp1(x,q_new,x_new)';
    plot(q_new,'k-')
    ylim([0 3])
    pause(0.1);
end

end