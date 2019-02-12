close all 
clear all

%%%%%%%%%%%%%%%%%%%%%%%% Q1 %%%%%%%%%%%%%%%%%%%%%%%%
%load dicom image
info = dicominfo('image18.dcm'); 
Y = dicomread(info); 
%convert to double
X = im2double(Y); 
%scale pixel values
I = (X/max(X(:)))*255; 
figure; imshow(I,[]);

%%%%%%%%%%%%%%%%%%%%%%%% Q2 %%%%%%%%%%%%%%%%%%%%%%%%
%crop 128x128 selection
crop_size = 128;
rect = [150 150 127 127]; 
cropped = imcrop(I,rect); 

%plot cropped image
figure; 
imshow(cropped,[]); 
title('Original Cropped Image');

%%%%%%%%%%%%%%%%%%%%%%%% Q3 %%%%%%%%%%%%%%%%%%%%%%%%
%pad for rotation
padder = ceil(length(cropped)*0.2);
I_padded = padarray(cropped,[padder padder],'both');
%rotate 0-179 degrees in 1 degree increments 
n_projections = 180; 
delta_theta = 1; 
%%get projections at each degree rotation
for i=1:n_projections
    pic = imrotate(I_padded,(i-1)*delta_theta,'crop');
    sinogram(:,i) = sum(pic); 
end

%plot sinogram
figure;
imshow(sinogram, []);
title('Sinogram');
xlabel('\theta (degrees)');
ylabel('X\prime');

%%%%%%%%%%%%%%%%%%%%%%%% Q4 %%%%%%%%%%%%%%%%%%%%%%%%
% 1D Fourier Transform of projections 
shifted = fftshift(sinogram); 
FT1 = fftshift(fft(shifted)); 
figure; 
imshow(log(abs(FT1)),[]);
title('1D FFT of Sinogram'); 

%%%%%%%%%%%%%%%%%%%%%%%% Q5 %%%%%%%%%%%%%%%%%%%%%%%%
%assemble 2-d transform from 1-d projections
angles = repmat((0:180-1)*pi/180,length(FT1),1); 
matrix_bound = (length(FT1)-1)/2; 
rhos = repmat((-ceil(matrix_bound):floor(matrix_bound))',1,n_projections); 
[x,y] = pol2cart(angles,rhos);
[XI,YI] = meshgrid((-ceil(matrix_bound):floor(matrix_bound))); 

%%%%%%%%%%%%%%%%%%%%%%%% Q6 %%%%%%%%%%%%%%%%%%%%%%%%
FT2_assembled = griddata(x,y,FT1,XI,YI,'linear'); 
FT2_assembled(isnan(FT2_assembled)) = 0; 

%plot constructed 2D fft 
figure; 
imshow(log(abs(FT2_assembled)),[]); 
title('Constructed 2D FFT from 1D FFT'); 

%%%%%%%%%%%%%%%%%%%%%%%% Q7 %%%%%%%%%%%%%%%%%%%%%%%%
%direct 2-d fourier transform of image 
FT2 = fftshift(fft2(fftshift(I_padded))); 


%%%%%%%%%%%%%%%%%%%%%%%% Q8 %%%%%%%%%%%%%%%%%%%%%%%%
%plot magnitude of 2D FFT
figure;
imshow(log(abs(FT2)),[]);  
title('Direct 2D FFT'); 

%%%%%%%%%%%%%%%%%%%%%%%% Q9 %%%%%%%%%%%%%%%%%%%%%%%%
%take inverse 2D FFT to plot reconstructed image 
inverse_fft2 = abs(fftshift(ifft2(fftshift(FT2)))); 
inverse_cropped = imcrop(inverse_fft2,[padder+1 padder+1 crop_size-1 crop_size-1]);
figure, imshow(inverse_cropped,[]); 
title('Reconstructed Image from Inverse 2D FFT');




