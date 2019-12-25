%written by: Jui-Wen, Yeh 
%For PePe's second lab
%try to analyze the picture, distinguish colors from it
clc; clear all
pic = imread ('SharedScreenshot.jpg'); 
figure (1)
imshow(pic)
gray_pic = rgb2gray(pic); 
[size_pic_x, size_pic_y] = size (gray_pic); 
R_shell_pic = uint8(pic(:, :, 1));
G_shell_pic = uint8(pic(:, :, 2));
B_shell_pic = uint8(pic(:, :, 3)); 

% test out substract method
R = imsubtract(R_shell_pic, gray_pic); 
figure (200)
R = im2bw(R, 0.32);
imshow(R)

%do the video processing 
video_1 = VideoReader('video-1577262358.mp4'); 
get (video_1) 

%cut video into pcitures and substract red dot
NumberOfFrames =uint8( video_1.Duration .* video_1.FrameRate); 
video_tear_down = zeros(video_1.Height, video_1.Width, NumberOfFrames);
video_tear_down_process = video_tear_down; 
for img = 1:NumberOfFrames
    process  = read(video_1, img); 
    %process to find the red dot in the image
    gray_process = rgb2gray(process);
    R_shell_process = uint8(process(:, :, 1));    
    R_process = imsubtract(R_shell_process, gray_process); 
    R_process = im2bw(R_process, 0.32);
    video_tear_down_process(:, :, img) = R_process; 
    figure(500)
    imshow(R_process)
end 
