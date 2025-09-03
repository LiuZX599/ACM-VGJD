close all;
clear all;

img_num = 2; % Input image number (1~5)

%% Read input image
Img = imread([num2str(img_num), '.bmp']);
Img = double(Img(:,:,1));
fmin  = min(Img(:));
fmax  = max(Img(:));
Img = 255*(Img-fmin)/(fmax-fmin);

params = image_params(img_num);

%% Median filtering
midfilter = medfilt2(Img, [5 5], 'symmetric');
Imgf = (1 - params.s) .* Img + params.s .* midfilter; % Eq.10

%% Initialize contour
c0 = 1;
initialLSF = ones(size(Imgf(:,:,1))).*c0;
initialLSF(params.lsf_row(1):params.lsf_row(2), params.lsf_col(1):params.lsf_col(2)) = -c0;
u = initialLSF;

%% Show initial contour
figure;imagesc(Img,[0, 255]);colormap(gray);hold on;axis off,axis equal
title('Initial contour');
[c,h] = contour(u,[0 0],'g','LineWidth',1);
pause(0.1);

%% Parameter settings
sigma = params.sigma; %Initial Gaussian kernel scale
iter_outer = params.iter_outer; % Number of iterations
iter_inner = 10;   
timestep = 1;
epsilon = 1;
k = params.k; % Mean filter window size
G = fspecial('average', k);
alfa = params.alfa; % Driving force 
min_sigma = 10;
beta = std2(Imgf); 

b=ones(size(Imgf));  

%% Level set function iteration
for n=1:iter_outer 
    if sigma>min_sigma
        if (n ~= 0 && mod(n, 2) == 0)
            sigma = (sigma - min_sigma) / 2 * (exp(-0.025 * n) + exp(-0.001 * n)) + min_sigma; % Eq.15
        end
    else
        sigma =sigma;
    end

    if(timestep < 1 && n~= 0)
        timestep = timestep + 0.01; % Eq.16
    end

    K=fspecial('gaussian',round(2*sigma)*2+1,sigma);
    [u, b, C]= VGJD(u,Imgf, b, K, timestep, epsilon, iter_inner,alfa,beta);
    u = (8*u)./sqrt(1+(8*u).^2); % Eq.17
    u = imfilter(u,G,'symmetric'); % Eq.18

    if mod(n,1)==0
        pause(0.001);
        imagesc(Img,[0, 255]); colormap(gray); axis off; axis equal;
        hold on;
        contour(u,[0 0],'r');
        iterNum=[num2str(n), ' iterations'];
        title(iterNum);
        hold off;
    end
   
end

%% Display final segmentation result
imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
[c,h] = contour(u,[0 0],'r','LineWidth',1);
totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);

%% Visualize final level set function
figure;
mesh(u);
colorbar;
title('Final level set function');


