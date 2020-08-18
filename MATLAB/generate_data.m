%% 1. Configuration
save_dir = 'D:/DATA/'; % directory where the data will be saved.
signal_flag = true;  % true for signal-present (g1), false for signal-absent (g0) images.
Nimg = 500; % number of 2D images (per noise type) to be generated and saved at a time.
Nset = 8; % number of image sets (each set has Nimg images). 
diameter = 2; % signal diameter
VGF = 50; % volume glandular fraction
ROIx = 129; ROIy = 129; ROIz = 129; tt = 64; % pixel size of a reconstructed volume 

warning('off','all')
%% 2. Filter
[X,Y,Z] = meshgrid(1:ROIx,1:ROIy,1:ROIz);
idx = sqrt((X-floor((ROIx)/2+1)).^2+(Y-floor((ROIy)/2+1)).^2+(Z-floor((ROIz)/2+1)).^2)<=tt;
W = single(idx);

ROI_filter = 257; 

kf = -floor(ROI_filter/2) : floor(ROI_filter/2);
lf = -floor(ROI_filter/2) : floor(ROI_filter/2); 
mf = -floor(ROI_filter/2) : floor(ROI_filter/2);

pixel_size = 0.388/2;
delta_f_filter = 1/(ROI_filter*pixel_size);

[X,Y,Z] = meshgrid(kf, lf, mf);
f_filter = delta_f_filter*sqrt(X.^2+Y.^2+Z.^2);

beta = 3;
filter = 1./(f_filter.^(beta/2));
filter(floor(ROI_filter/2)+1,floor(ROI_filter/2)+1,floor(ROI_filter/2)+1) = 0;
filter(floor(ROI_filter/2)+1,floor(ROI_filter/2)+1,floor(ROI_filter/2)+1) = 2*max(max(max(filter)));

%% 3. Image Generation
hann_linear_trans  = zeros(ROIx,ROIx,Nimg);
shepp_linear_trans = zeros(ROIx,ROIx,Nimg);
ramp_linear_trans  = zeros(ROIx,ROIx,Nimg);
hann_linear_longi  = zeros(ROIx,ROIx,Nimg);
shepp_linear_longi = zeros(ROIx,ROIx,Nimg);
ramp_linear_longi  = zeros(ROIx,ROIx,Nimg);

hann_sinc_trans  = zeros(ROIx,ROIx,Nimg);
shepp_sinc_trans = zeros(ROIx,ROIx,Nimg);
ramp_sinc_trans  = zeros(ROIx,ROIx,Nimg);
hann_sinc_longi  = zeros(ROIx,ROIx,Nimg);
shepp_sinc_longi = zeros(ROIx,ROIx,Nimg);
ramp_sinc_longi  = zeros(ROIx,ROIx,Nimg);


for jj=1:Nset
    for ii=1:Nimg
        tic 

        % 1. Breast Volume
        bg_temp = ifftn(ifftshift(fftshift(fftn(randn(ROI_filter,ROI_filter,ROI_filter))).*filter));
        img_temp = bg_temp(floor(ROI_filter/2)+1-tt:floor(ROI_filter/2)+1+tt,floor(ROI_filter/2)+1-tt:floor(ROI_filter/2)+1+tt,floor(ROI_filter/2)+1-tt:floor(ROI_filter/2)+1+tt);
        temp = img_temp;
        sorted = sort(temp(:),'descend');
        threshold_val = sorted(round(length(sorted)*(VGF/100)));

        idx = temp > threshold_val;
        img = single(idx)*0.233 + single(~idx)*0.194;

        if signal_flag == true
            % Signal insertion
            [X,Y,Z] = meshgrid(1:ROIx,1:ROIy,1:ROIz);
            idx = sqrt((X-floor(ROIx/2+1)).^2+(Y-floor(ROIy/2+1)).^2+(Z-floor(ROIz/2+1)).^2) < diameter/pixel_size/2;
            img(idx) = 0.238;

            img = img.*W;
        end

        % 2. Projection
        ParamSetting;
        proj = zeros(param.nu, param.nv, param.nProj,'single');
        parfor i = 1:param.nProj    
            proj(:,:,i) = projection(img,param,i);
        end

        ParamSetting;

        N0=6914; %
        Nin=N0*exp(zeros(size(proj)));
        proj=proj+10*(log(Nin)-log(poissrnd(N0*ones(size(proj))))); % x10 because our unit is mm but unit of mu is cm


        % 3. Filtered Projection
        param.filter='hann'; 
        proj_filtered_hann  = filtering(proj,param); % vertical: y, horizontal: z
        param.filter='shepp-logan';
        proj_filtered_shepp = filtering(proj,param);
        param.filter='ram-lak';
        proj_filtered_ramp  = filtering(proj,param);


        % 4. Sinc Interpolation
        IF = 8;
        center = floor((IF*param.nu-(IF-1))/2+1);
        L = length(center-floor(IF*param.nu/2.8):center+floor(IF*param.nu/2.8));
        proj_filtered_hann_sinc3  = zeros(L,L,param.nProj);
        proj_filtered_shepp_sinc3  = zeros(L,L,param.nProj);
        proj_filtered_ramp_sinc3  = zeros(L,L,param.nProj);

        parfor iview = 1:param.nProj
            proj_filtered_hann_sinc1 = interpft(proj_filtered_hann(:,:,iview),param.nu*IF);
            proj_filtered_hann_sinc2 = interpft(proj_filtered_hann_sinc1,param.nu*IF,2);
            proj_filtered_hann_sinc3(:,:,iview) = proj_filtered_hann_sinc2(center-floor(IF*param.nu/2.8):center+floor(IF*param.nu/2.8),center-floor(IF*param.nu/2.8):center+floor(IF*param.nu/2.8));

            proj_filtered_shepp_sinc1 = interpft(proj_filtered_shepp(:,:,iview),param.nu*IF);
            proj_filtered_shepp_sinc2 = interpft(proj_filtered_shepp_sinc1,param.nu*IF,2);
            proj_filtered_shepp_sinc3(:,:,iview) = proj_filtered_shepp_sinc2(center-floor(IF*param.nu/2.8):center+floor(IF*param.nu/2.8),center-floor(IF*param.nu/2.8):center+floor(IF*param.nu/2.8));   

            proj_filtered_ramp_sinc1 = interpft(proj_filtered_ramp(:,:,iview),param.nu*IF);
            proj_filtered_ramp_sinc2 = interpft(proj_filtered_ramp_sinc1,param.nu*IF,2);
            proj_filtered_ramp_sinc3(:,:,iview) = proj_filtered_ramp_sinc2(center-floor(IF*param.nu/2.8):center+floor(IF*param.nu/2.8),center-floor(IF*param.nu/2.8):center+floor(IF*param.nu/2.8));        
        end

        % 5. Reconstruction
        ParamSetting1;
        tempx = (1:1/IF:129); tempy = (1:1/IF:129);
        Recon_hann_linear = 0; Recon_shepp_linear = 0; Recon_ramp_linear = 0;  
        Recon_hann_sinc = 0;   Recon_shepp_sinc = 0;   Recon_ramp_sinc = 0;
        [x,y] = meshgrid(tempx(center-floor(IF*param.nu/2.8):center+floor(IF*param.nu/2.8)),tempy(center-floor(IF*param.nu/2.8):center+floor(IF*param.nu/2.8)));
        [xx,yy] = meshgrid(param.xs,param.ys);

        parfor iview = 1:param.nProj  
            angle_rad = param.deg(iview)/360*2*pi;


            rx = xx.*cos(angle_rad-pi/2) + yy.*sin(angle_rad-pi/2);
            ry = -xx.*sin(angle_rad-pi/2) + yy.*cos(angle_rad-pi/2);

            pu = single(((rx.*(param.DSD)./(ry + param.DSO))+param.us(1))/(-param.du) + 1);
            Ratio = (single(param.DSO.^2./(param.DSO+ry).^2));    

            vol_hann_linear  = zeros(param.nx,param.ny,param.nz,'single');
            vol_shepp_linear = zeros(param.nx,param.ny,param.nz,'single');
            vol_ramp_linear  = zeros(param.nx,param.ny,param.nz,'single');
            vol_hann_sinc    = zeros(param.nx,param.ny,param.nz,'single');
            vol_shepp_sinc   = zeros(param.nx,param.ny,param.nz,'single');
            vol_ramp_sinc    = zeros(param.nx,param.ny,param.nz,'single');

            for iz = 1:param.nz   
                pv = single(((param.zs(iz)*(param.DSD)./(ry + param.DSO))-param.vs(1))/param.dv+1);
                vol_hann_linear(:,:,iz)  = (Ratio.*interp2(proj_filtered_hann(:,:,iview)',pu,pv,param.interptype));
                vol_shepp_linear(:,:,iz) = (Ratio.*interp2(proj_filtered_shepp(:,:,iview)',pu,pv,param.interptype));
                vol_ramp_linear(:,:,iz)  = (Ratio.*interp2(proj_filtered_ramp(:,:,iview)',pu,pv,param.interptype));
                vol_hann_sinc(:,:,iz)    = (Ratio.*interp2(x,y,proj_filtered_hann_sinc3(:,:,iview)',pu,pv,param.interptype));
                vol_shepp_sinc(:,:,iz)   = (Ratio.*interp2(x,y,proj_filtered_shepp_sinc3(:,:,iview)',pu,pv,param.interptype));
                vol_ramp_sinc(:,:,iz)    = (Ratio.*interp2(x,y,proj_filtered_ramp_sinc3(:,:,iview)',pu,pv,param.interptype));
            end

            vol_hann_linear(isnan(vol_hann_linear))=0;
            vol_shepp_linear(isnan(vol_shepp_linear))=0;
            vol_ramp_linear(isnan(vol_ramp_linear))=0;
            vol_hann_sinc(isnan(vol_hann_sinc))=0;
            vol_shepp_sinc(isnan(vol_shepp_sinc))=0;
            vol_ramp_sinc(isnan(vol_ramp_sinc))=0;

            Recon_hann_linear = Recon_hann_linear + vol_hann_linear;
            Recon_shepp_linear = Recon_shepp_linear + vol_shepp_linear;
            Recon_ramp_linear = Recon_ramp_linear + vol_ramp_linear;
            Recon_hann_sinc = Recon_hann_sinc + vol_hann_sinc;
            Recon_shepp_sinc = Recon_shepp_sinc + vol_shepp_sinc;
            Recon_ramp_sinc = Recon_ramp_sinc + vol_ramp_sinc;

        end

        hann_linear_trans(:,:,ii)  = Recon_hann_linear(:,:,65);
        shepp_linear_trans(:,:,ii) = Recon_shepp_linear(:,:,65);
        ramp_linear_trans(:,:,ii)  = Recon_ramp_linear(:,:,65);       
        hann_linear_longi(:,:,ii)  = squeeze(Recon_hann_linear(65,:,:));
        shepp_linear_longi(:,:,ii) = squeeze(Recon_shepp_linear(65,:,:));
        ramp_linear_longi(:,:,ii)  = squeeze(Recon_ramp_linear(65,:,:));

        hann_sinc_trans(:,:,ii)  = Recon_hann_sinc(:,:,65);
        shepp_sinc_trans(:,:,ii) = Recon_shepp_sinc(:,:,65);
        ramp_sinc_trans(:,:,ii)  = Recon_ramp_sinc(:,:,65);       
        hann_sinc_longi(:,:,ii)  = squeeze(Recon_hann_sinc(65,:,:));
        shepp_sinc_longi(:,:,ii) = squeeze(Recon_shepp_sinc(65,:,:));
        ramp_sinc_longi(:,:,ii)  = squeeze(Recon_ramp_sinc(65,:,:));

        clc

        ii
        toc
    end

    save([save_dir,'g',num2str(signal_flag),'_',num2str(jj)], ...
             'hann_linear_trans','shepp_linear_trans','ramp_linear_trans', ...
             'hann_linear_longi','shepp_linear_longi','ramp_linear_longi', ...
             'hann_sinc_trans','shepp_sinc_trans','ramp_sinc_trans', ...
             'hann_sinc_longi','shepp_sinc_longi','ramp_sinc_longi');
end

%% 4. Plotting Sample Images
close all
dspmin = 0.12; dspmax = 0.30;
figure;
ii = 1;
subplot(261); imshow(hann_linear_trans(:,:,ii),[]);
subplot(262); imshow(shepp_linear_trans(:,:,ii),[]);
subplot(263); imshow(ramp_linear_trans(:,:,ii),[]);
subplot(264); imshow(hann_sinc_trans(:,:,ii),[]);
subplot(265); imshow(shepp_sinc_trans(:,:,ii),[]);
subplot(266); imshow(ramp_sinc_trans(:,:,ii),[]);
subplot(267); imshow(hann_linear_longi(:,:,ii)',[]);
subplot(268); imshow(shepp_linear_longi(:,:,ii)',[]);
subplot(269); imshow(ramp_linear_longi(:,:,ii)',[]);
subplot(2,6,10); imshow(hann_sinc_longi(:,:,ii)',[]);
subplot(2,6,11); imshow(shepp_sinc_longi(:,:,ii)',[]);
subplot(2,6,12); imshow(ramp_sinc_longi(:,:,ii)',[]);

%% 5. Formatting Datasets
N_names = {'hann_linear_trans', 'shepp_linear_trans',   'ramp_linear_trans', ...
           'hann_sinc_trans',   'shepp_sinc_trans',     'ramp_sinc_trans', ...
           'hann_linear_longi', 'shepp_linear_longi',   'ramp_linear_longi',...
           'hann_sinc_longi',   'shepp_sinc_longi',     'ramp_sinc_longi'};

files = dir([save_dir,'g0_*.mat']);

for N = 1:12
    g0 = [];
    for i = 1:length(files)
        load([save_dir,files(i).name], N_names{N});
        switch N
            case 1
                g_temp = hann_linear_trans;
            case 2
                g_temp = shepp_linear_trans;
            case 3
                g_temp = ramp_linear_trans;
            case 4
                g_temp = hann_sinc_trans;            
            case 5
                g_temp = shepp_sinc_trans;
            case 6
                g_temp = ramp_sinc_trans;
            case 7
                g_temp = hann_linear_longi;
            case 8
                g_temp = shepp_linear_longi;                  
            case 9
                g_temp = ramp_linear_longi;
            case 10
                g_temp = hann_sinc_longi;
            case 11
                g_temp = shepp_sinc_longi;
            case 12
                g_temp = ramp_sinc_longi; 
        end
        g0 = cat(3,g0,g_temp);
    end
    save([save_dir,'N',num2str(N),'_g0.mat'], 'g0', '-v7.3');
end

files = dir('D:/0902/g1_*.mat');

for N = 1:12
    g1 = [];
    for i = 1:length(files)
        load([save_dir,files(i).name], N_names{N});
        switch N
            case 1
                g_temp = hann_linear_trans;
            case 2
                g_temp = shepp_linear_trans;
            case 3
                g_temp = ramp_linear_trans;
            case 4
                g_temp = hann_sinc_trans;            
            case 5
                g_temp = shepp_sinc_trans;
            case 6
                g_temp = ramp_sinc_trans;
            case 7
                g_temp = hann_linear_longi;
            case 8
                g_temp = shepp_linear_longi;                  
            case 9
                g_temp = ramp_linear_longi;
            case 10
                g_temp = hann_sinc_longi;
            case 11
                g_temp = shepp_sinc_longi;
            case 12
                g_temp = ramp_sinc_longi; 
        end
        g1 = cat(3,g1,g_temp);
    end
    save([save_dir,'N',num2str(N),'_g1.mat'], 'g1', '-v7.3');
end
