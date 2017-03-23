function [ output ] = get_inlier( error , t )

output=[];
for i=1:size(error,1)
    if(error(i)<=t)
        output = [output,i];
    end
end;

end

