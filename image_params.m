function params = image_params(img_num)
    switch img_num
        case 1  % ICT image with no obvious features
            params.s = 0;
            params.lsf_row = [90, 240];   
            params.lsf_col = [140, 270];  
            params.sigma = 40;
            params.iter_outer = 30;
            params.k = 7;
            params.alfa = 3;

        case 2  % ICT image with scattering artifacts
            params.s = 0;
            params.lsf_row = [221, 294];  
            params.lsf_col = [230, 315]; 
            params.sigma = 60;
            params.iter_outer = 35;
            params.k = 7;
            params.alfa = 3;
            
        case 3 % ICT image with severe noise
            params.s = 0.4;
            params.lsf_row = [94, 600];  
            params.lsf_col = [99, 707];  
            params.sigma = 70;
            params.iter_outer = 50;
            params.k = 7;
            params.alfa = 1.5;

        case 4  % Initial contour deviates from target
            params.s = 0;
            params.lsf_row = [242, 312];   
            params.lsf_col = [30, 100];  
            params.sigma = 80;
            params.iter_outer = 40;
            params.k = 9;
            params.alfa = 1;

        case 5  % ICT image with severe scattering artifacts and noise
            params.s = 0.3;
            params.lsf_row = [83, 579];   
            params.lsf_col = [112, 415];  
            params.sigma = 105;
            params.iter_outer = 100;
            params.k = 9;
            params.alfa = 0.4;
            
    end
end