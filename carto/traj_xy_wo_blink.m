function traj_xy_wo_blink(tab_param, color)

global N_PARAM PARAM_ALPHA PARAM_I PARAM_J

%% --- go through traces ---
Ntrc = size(tab_param, 2);

disp('traj :            ')

for itrc = 1:Ntrc
    alpha = tab_param((PARAM_ALPHA-1):N_PARAM:end, itrc);
    t = find(alpha>0);
    pos_y = tab_param(N_PARAM*(t-1)+PARAM_I-1, itrc); % !! i=y & x=j !!
    pos_x = tab_param(N_PARAM*(t-1)+PARAM_J-1, itrc);
    plot(pos_x+1, pos_y+1, color)
    
    if mod(itrc,500)==0
        fprintf([repmat('\b', 1, 11) '%5i/%5i'], itrc, Ntrc)
        drawnow
    end
end

fprintf([repmat('\b', 1, 11) '%5i/%5i\n'], Ntrc, Ntrc)

%%%