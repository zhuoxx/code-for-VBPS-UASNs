function [vut] = valueofinformation(vusenddelay_t,betau,vu0,alphau,Tu)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

%% vu
if vusenddelay_t<=Tu
   vut=betau*vu0+(1-betau)*vu0*exp(-(vusenddelay_t)/alphau);
else
    vut=0;
end


end

