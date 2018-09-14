function q_new = fe_adv_sla(delta_t,u_func,x,q_init,target_time)

q_new = q_init;

for i = 1:(target_time/delta_t)
    u = u_func(i*delta_t);
    x_new = x + -u*delta_t; % new set of x values
    id = x_new < 0;
    x_new(id) = 1 - abs(x_new(id));
    id = x_new > 1;    
    x_new(id) = x_new(id)-1;
%     q_new = interp1([flipud(-x(2:end));x;1+x],[flipud(q_new(2:end));q_new;q_new],q_new,x_new)';
    q_new = interp1([-x(2:end);x;1+x],[q_new(2:end);q_new;q_new],x_new,'spline');
    plot(x,q_new,'k-')
    ylim([0 3])
    pause(0.1);
end

end