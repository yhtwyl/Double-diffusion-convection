function [ d ] = func_vel1( f_name,var )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
  if (var == 'u')
    end_x = 100;end_y = 100;end_z = 100;
  elseif(var == 'v')
    end_x = 100;end_y = 100;end_z = 100;
  elseif(var == 'w')
    end_x = 100;end_y = 100;end_z = 100;
  else
    end_x = 100;end_y = 100;end_z = 100;
  end
fid = fopen(f_name,'r');
data = fread(fid,'real*4');
size(data);
%    d = zeros(101,101,101);
      for k = 1:end_z                               
        for j = 1:end_y
            for i = 1:end_x
                d(i,j,k) = data(i + (j-1)*100 + 100*100*(k-1));
            end
        end
%%%	    k,d(:,:,k)
      end
end
