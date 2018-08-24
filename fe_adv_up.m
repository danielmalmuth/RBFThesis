function q_new = fe_adv_up(delta_t,u,q_init,target_time,h)

q_new = q_init;

for i = 1:(target_time/delta_t)
    Dq = up_adv_diff(u,q_new,h);
    q_new = q_new + (delta_t*u).*Dq;
%     plot(q_new,'k-')
%     ylim([0 1])
%     pause(0.2);
end

end