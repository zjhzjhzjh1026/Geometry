function [ output ] = get_inlier( error , t )

output=[];
for i=1:60
    if(error(i)<=t)
        output = [output,i];
    end
end;

end

