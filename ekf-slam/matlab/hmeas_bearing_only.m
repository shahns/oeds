function [zpred, Hjacobian] = hmeas_bearing_only(xpred, state_ind)

%% initializations
p = xpred(1:2);
no_visible_lm = length(state_ind);
zpred = zeros(no_visible_lm, 1);
Hjacobian = zeros(no_visible_lm, length(xpred));
%% predict
for lm = 1:no_visible_lm
    xm = xpred(2*state_ind(lm)-1);
    ym = xpred(2*state_ind(lm));
    zpred(lm) = atan2d(ym-p(2), xm-p(1));
    Hjacobian(lm, 1) = -(ym - p(2)) / ((xm - p(1))^2 + (ym - p(2))^2);
    Hjacobian(lm, 2) = (xm - p(1)) / ((xm - p(1))^2 + (ym - p(2))^2);
    Hjacobian(lm, 2*state_ind(lm)-1) = (ym - p(2)) / ((xm - p(1))^2 + (ym - p(2))^2);
    Hjacobian(lm, 2*state_ind(lm)) = -(xm - p(1)) / ((xm - p(1))^2 + (ym - p(2))^2);
    
end
    Hjacobian = (pi / 180) * Hjacobian;
end
