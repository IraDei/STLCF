function [res] = imGrayNorm(imin, pos_only)
%IMGRAYNORM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
gmin = min(imin,[],'all');
gmax = max(imin,[],'all');
if ~pos_only
    gdel = gmax - gmin + .000001;
    res = (imin - gmin)./gdel;
    fprintf('input image: min/max = %g/%g\n', gmin, gmax);
else
    res = exp(imin)./exp(gmax);
end
end

