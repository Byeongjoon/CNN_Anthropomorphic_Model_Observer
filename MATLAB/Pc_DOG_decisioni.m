function [tt] = Pc_DOG_decisioni(g1,g0,ROI,cutoff,alpha,Q,Nc,Nimg_train,g_test,int_noise)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    c = zeros(ROI,ROI,Nc); t = zeros(ROI*ROI,Nc);
    for ch=1:Nc
        c(:,:,ch) = DOG(ROI,cutoff,alpha,ch,Q);
        t_temp = fftshift(ifftn(ifftshift(c(:,:,ch))));
        t(:,ch) = t_temp(:);
    end
    T = t'; 

    % Training
    X0 = zeros(Nc, Nimg_train); X1 = zeros(Nc, Nimg_train);
    for i=1:Nimg_train
        temp1 = g1(:,:,i);   temp0 = g0(:,:,i);
        temp11 = T*temp1(:); temp00 = T*temp0(:);
        X1(:,i) = temp11;     X0(:,i) = temp00;
    end
    K1 = cov(X1',1); K0 = cov(X0',1);
    K = 0.5*(K1+K0);

    inv_K = inv(K);


    delta_s = mean(g1(:,:,1:Nimg_train),3)-mean(g0(:,:,1:Nimg_train),3);
    delta_v = T*delta_s(:);
    w = delta_v'*(inv_K);
    
    % test

    for i=1:length(g_test(1,1,:))
        g_test_temp = g_test(:,:,i);
        tt(i) = w*(T*g_test_temp(:))+int_noise*randn(1);
    end
end
    


