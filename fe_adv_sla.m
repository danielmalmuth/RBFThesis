function q_new = fe_adv_sla(delta_t,u_func,x,q_init,target_time)

q_new = q_init;

for i = 1:round((target_time/delta_t))
    u = u_func(i*delta_t,x);
    x_new = x + (delta_t*-u); % new set of x values
    id = x_new > 1;
    x_new(id) = x_new(id)-1;
    id = x_new < 0;
    x_new(id) = x_new(id)+1;
    q_new = interp1([-1+x(2:end);x;1+x(2:end)],...
        [q_new(2:end);q_new;q_new(2:end)],x_new,'spline');
%     plot(x,q_new,'k-')
%     ylim([0 3])
%     pause(0.1);
end

end