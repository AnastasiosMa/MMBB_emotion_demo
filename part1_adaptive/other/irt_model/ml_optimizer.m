function [th, th_idx] = ml_optimizer(th,type,u,j,p_correct,p_incorrect)
optimizer = @(th_idx,u,j) sum(u-p_correct(j,th_idx)')/(-sum(p_correct(j,th_idx)'.*p_incorrect(j,th_idx)'));
init_th = 0; %initial starting point of theta
min_criterion = 0.05; %optimizer stops if Δθ < criterion
n_iter = 500; %number of iterations before manual stop
theta_low = -6;
theta_high = 6;
theta_step = 0.05;
theta_range = round(theta_low:theta_step:theta_high,2);
delta_th = 1;
diffulty_change = 1;
iter = 1;
if type==2
    th = 0;
    while abs(delta_th)>min_criterion
        if iter>n_iter
            break
        end
        if th<theta_low
            th_idx = 1;
            th = theta_low;
            if ~any(u) %if all answers wrong, break
                break
            end
        elseif th>theta_high
            th_idx = length(theta_range);
            th = theta_high;
            if all(u)%if all answers correct, break
                break
            end
        else
            %find discrete index for theta within range
            th_idx = find(th<=theta_range+theta_step & th>theta_range);
        end
        if isempty(th_idx)
            keyboard
        end
        delta_th = optimizer(th_idx,u,j);
        th = th-delta_th;
        iter = iter+1;
    end
end
if type==1
   if u(end)==1
      th = th + (theta_high-th)/(2*theta_high/diffulty_change);
   else
      th = th + (theta_low-th)/(2*theta_high/diffulty_change);
   end
end
if th<theta_low
    th_idx = 1;
elseif th>theta_high
    th_idx = length(theta_range);
else
    th_idx = find(th<=theta_range+theta_step & th>theta_range);
end
if isnan(th) | isempty(th_idx)
                keyboard
            end
end