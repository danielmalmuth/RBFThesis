function q_new = rk4_adv_sla(delta_t,u_func,x,q_init,target_time)

q_new = q_init;
y = zeros(length(x),1);

for i = 1:round((target_time/delta_t))
    k1 = -delta_t.*u_func(i*delta_t,x);
    k2 = -delta_t.*u_func(i*delta_t + -(delta_t/2),x + (k1/2));
    k3 = -delta_t.*u_func(i*delta_t + -(delta_t/2),x + (k2/2));
    k4 = -delta_t.*u_func(i*delta_t + -(delta_t),x + k3);
    x_new = x + (k1 + 2.*k2 + 2.*k3 + k4)./6;
    id = x_new > 1;
    x_new(id) = x_new(id)-1;
    id = x_new < 0;
    x_new(id) = x_new(id)+1;
    q_new = interp1([-1+x(2:end);x;1+x(2:end)],...
        [q_new(2:end);q_new;q_new(2:end)],x_new,'spline');
    
%     lam = rbffit_test(x,y,q_new,"Backslash",true,0,1);
%     q_new = rbfval2(lam,x,y,x_new,y,1);
%     plot(x,q_new,'k-')
%     ylim([0 3])
%     pause(0.1);
end

end