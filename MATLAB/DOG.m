function [c] = DOG(M,cutoff,alpha,j,Q)
% low pass noise structure

c = zeros(M,M);
law_k = zeros(M,M);
for x = 1 : M
    for y = 1 : M
        law_k(x,y) = (sqrt((floor(M/2+1)-x)^2+(floor(M/2+1)-y)^2))/(M/2);
        c(x,y)=exp(-0.5*(law_k(x,y)/(Q*cutoff*alpha^j))^2)-exp(-0.5*(law_k(x,y)/(cutoff*alpha^j))^2);
    end
end

end
% pixel_size = 0.1;
% X=linspace(-1/(pixel_size*2),1/(pixel_size*2),M); Y=linspace(-1/(pixel_size*2),1/(pixel_size*2),M);
% for x = 1 : M
%     for y = 1 : M
%         law_k(x,y) = sqrt(X(x)^2+Y(y)^2);
%         c(x,y)=exp(-0.5*(law_k(x,y)/(Q*cutoff*alpha^j))^2)-exp(-0.5*(law_k(x,y)/(cutoff*alpha^j))^2);
%     end
% end
