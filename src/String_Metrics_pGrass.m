function [StringOut] = String_Metrics_pGrass(StringIn, CaseName, Val1)
%  
% 

curr_string = sprintf('%30s %10.4f \n',CaseName,Val1);
StringOut   = strcat(StringIn,curr_string,'\n');

end
