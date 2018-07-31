function q_new = fe_adv(delta_t,u,D,q_init,target_time)

q_new = q_init;

for i = 1:(1/(delta_t*target_time))
    q_new = q_new + (delta_t*u).*D*q_new;
    plot(q_new,'k-')
    ylim([0 1])
    pause(0.2);
end

end