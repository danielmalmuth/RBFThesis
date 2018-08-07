function q_new = be_adv(delta_t,u,D,q_init,target_time)

q_new = q_init;

% Fix this
% L = (speye(size(D,1))-spdiags((delta_t*u))*D);
L = (speye(size(D,1))-(delta_t*u)*D);
decompL = decomposition(L);
for i = 1:(target_time/delta_t)
    q_new = decompL\q_new;
%     plot(q_new,'k-')
%     ylim([0 1])
%     pause(0.2);
end

end