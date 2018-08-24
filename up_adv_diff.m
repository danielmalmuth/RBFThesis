function Dq = up_adv_diff(u,q_new,h)

Dq = q_new;
n = length(q_new);

if u > 0
    Dq(1) = (q_new(1)-q_new(n))/h;
    Dq(n) = (q_new(n)-q_new(n-1))/h;
else
    Dq(1) = (q_new(2)-q_new(1))/h;
    Dq(n) = (q_new(1)-q_new(n))/h;
end

for i = 2:(n-1)
    if u > 0
        Dq(i) = (q_new(i)-q_new(i-1))/h;
    else
        Dq(i) = (q_new(i+1)-q_new(i))/h;
    end
end

end