function q_new = ...
    rk4_sla(delta_t,u_func,x,int_leng,q_init,target_time, ...
    rbf_exp,stencil,deg)

q_new = q_init;
q_new_ext = [q_new;q_new;q_new];

x_ext = [x-int_leng;x;x+int_leng];
tree = KDTreeSearcher(x_ext);

for i = 1:length(x)
    idx = knnsearch(tree,x(i),'k',stencil);
    S{i}.neighbors = idx;
    S{i}.d = localinterpmat_1D(x_ext(idx),rbf_exp,deg);
end
% map = repmat((1:length(x))',[3 1]);

for i = 1:round((target_time/delta_t))
    k1 = -delta_t.*u_func(i*delta_t,x);
    k2 = -delta_t.*u_func(i*delta_t + -(delta_t/2),x + (k1/2));
    k3 = -delta_t.*u_func(i*delta_t + -(delta_t/2),x + (k2/2));
    k4 = -delta_t.*u_func(i*delta_t + -(delta_t),x + k3);
    x_new = x + (k1 + 2*k2 + 2*k3 + k4)./6;
 
    for j = 1:length(x)
        idx = knnsearch(tree,x_new(j),'k',1);
        idx = mod(idx-1,length(x))+1;
        lam = S{idx}.d\[q_new_ext(S{idx}.neighbors);zeros(deg+1,1)];
        
%         lam = rbffit_test_1D(x_ext(idx),q_new_ext(idx),"Backslash",true,0,2); 
        if x_new(j) > (x_ext(2*length(x)) + x_ext(2*length(x)+1))/2
            x_new(j) = x_new(j) - x_ext(2*length(x)+1);
        elseif x_new(j) < (x_ext(length(x)+1) + x_ext(length(x)))/2
            x_new(j) = x_ext(2*length(x)+1) + x_new(j);
        end
        q_new(j) = ...
            localrbfinterp_1D(lam,x_new(j),x_ext(S{idx}.neighbors),...
            rbf_exp,deg); % Local method
        
%         q_new = rbfval2_1D(lam,x(idx),x_new,1); % Global method

    end
    q_new_ext = [q_new;q_new;q_new];

%     plot(x,q_new,'k-')
%     ylim([0 3])
%     pause(0.005);
end

end