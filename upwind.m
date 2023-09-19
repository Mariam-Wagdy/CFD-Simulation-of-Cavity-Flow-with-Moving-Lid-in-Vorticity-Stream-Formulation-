function [ex, ey]=upwind(u,v)
    %% upwind condition
    if u>0
        ex=1;
    elseif u<0
        ex=-1;
    else
        ex=0;
    end
    
    if v>0
        ey=1;
    elseif v<0
        ey=-1;
    else
        ey=0;
    end
end
