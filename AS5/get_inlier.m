function [ inlier , outlier ] = get_inlier( error , t )

inlier = [];
outlier = [];
for i=1:size(error,1)
    if(error(i)<=t)
        inlier = [inlier,i];
    else
        outlier = [outlier,i];
    end
end;

end

