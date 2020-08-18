function [t,w] = Pc_NPWE(g1,g0,ROI,eta,c,g_test)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


    law_k = zeros(ROI,ROI); 
    for x = 1 : ROI
        for y = 1 : ROI
            law_k(x,y) = (sqrt((floor(ROI/2+1)-x)^2+(floor(ROI/2+1)-y)^2))/(ROI)/0.1014;
        end
    end
    E = (law_k.^eta).*exp(-c*law_k.^2);
    delta_s = mean(g1,3)-mean(g0,3);
    
    w = ifftn(ifftshift((E.^2).*fftshift(fftn(delta_s))));

    % Test
     for i=1:length(g_test(1,1,:))
         g_test_temp = g_test(:,:,i);
         t(i) = w(:)'*g_test_temp(:);
    end
  
 

end
