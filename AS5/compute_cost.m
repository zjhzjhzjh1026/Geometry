function [ cost , inlier ] = compute_cost( error , t )

cost=0;
inlier = 0;
for i=1:size(error,1)
    if(error(i)<=t)
        cost = cost + error(i);
        inlier = inlier + 1;
    else
        cost = cost + t; 
    end
end;

end

