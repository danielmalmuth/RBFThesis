function q_new = be_adv(delta_t,u,D,q_init,target_time)

q_new = q_init;

for i = 1:(1/(delta_t*target_time))
    q_new = (eye(size(D,1))-(delta_t*u).*D)\q_new;
end

end