% Define colors in RGB space
black = [0, 0, 0];
red = [255, 0, 0];
cyan = [0, 255, 255];
green = [0,255,0];
blue = [0,0,255];
yellow = [255, 255, 0];
white = [255, 255, 255];
magenta =[255,0,255];

colors = {black,blue,cyan,green,red,yellow, magenta,white};

% Read text file
fileID = fopen('BarcodeText.txt', 'r');
formatspec = '%c';
A = fscanf(fileID, formatspec, Inf);
B = double(A);
C = dec2bin(B,9);
C = C - '0';
C = C';
BinaryInput = reshape(C, 1, []);

BitsPerChar = 3;
AmountOfBits = size(BinaryInput, 2);
RemainingBits = mod(AmountOfBits, 900);
FillBarcode = mod(RemainingBits, 30);





if RemainingBits ~= 0
    Addzeroes = 900 - RemainingBits;
    BinaryInput(1, AmountOfBits + Addzeroes) = 0;
else
end

% Create a directory for frames if it doesn't exist
if ~exist('frames', 'dir')
    mkdir('frames');
end


VideoUpload = VideoWriter('Test.avi','Uncompressed AVI');
VideoUpload.FrameRate = 2;  % Adjust the frame rate to a valid value

open(VideoUpload);

nFrames = ceil(AmountOfBits / (30 * 30 * BitsPerChar));
Header_maxframes = dec2bin(nFrames,9);
Header_maxframes=Header_maxframes-'0';

for frame = 1:nFrames
     BarcodeFrame = ones(30) * 9;

    for row = 1:30
        for col = 1:30
            bitIndex = (frame - 1) * 30 * 30 * BitsPerChar + (row - 1) * 30 * BitsPerChar + (col - 1) * BitsPerChar + 1;

            if bitIndex + BitsPerChar - 1 <= AmountOfBits
                colorBits = BinaryInput(bitIndex:bitIndex + BitsPerChar - 1);
                colorIdx = bi2de(colorBits, 'left-msb') + 1;

            else
                colorIdx = 1; % default to black if not enough bits
            end

            BarcodeFrame(row, col) = colorIdx;
        end
    end

% Create colored frame
ColoredFrame = zeros(30, 30, 3);
for row = 1:30
    for col = 1:30
        ColoredFrame(row, col, :) = colors{BarcodeFrame(row, col)};
    end
end

% Add border
WhiteSpace = padarray(ColoredFrame, [1 1 0], 1, 'both');
BarcodeBorder = padarray(WhiteSpace, [2 2 0], 'both');

% modify

Error_rem = rem(frame,2);
if Error_rem==0
Errorbits = [1 0 0];
else
Errorbits = [0 1 0];
end
Frame_num = dec2bin(frame,9);
Frame_num = Frame_num-'0';
Headerinfo = [Errorbits Frame_num Header_maxframes];
HeaderAmountOfBits = size (Headerinfo,2);
Header = zeros(1,30);
for row = 1:1
    for col = 1:30
    if col >= 23 && col <= 30
        % Display colors from the 'colors' array in the last 8 columns
        colorIdx2 = col - 22;
    else
        index = (col - 1) * BitsPerChar + 1;
        if index + BitsPerChar - 1 <= HeaderAmountOfBits
            colorBits2 = Headerinfo(index:index + BitsPerChar - 1);
            colorIdx2 = bi2de(colorBits2, 'left-msb')+1 ;
        else
         colorIdx2 = 1; % default to black for future use  
        end
    end
    Header(row, col) = colorIdx2;
    end
end

%Create colored header
ColoredHeader = zeros(1, 30, 3);
for row = 1:1
    for col = 1:30
        ColoredHeader(row, col, :) = colors{Header(row, col)};
    end
end

PadHdr = padarray(ColoredHeader,[1 1 0],1,'both');
BlackPadd = padarray(PadHdr, [0 2 0], 'both');
BlackPadd = padarray(BlackPadd, [2 0 0], 'post');
% Resize ColoredHeader to match the size of BarcodeBorder
%%desiredWidth = 42;
%%ColoredHeaderResized = imresize(PadHdr, [size(PadHdr, 1), desiredWidth]);


Hdr_barcode = cat(1, BarcodeBorder, BlackPadd);


PadBarcode = padarray(Hdr_barcode, [15 15 0], 1, 'both');
% Expand each pixel
blockSize = 15; % Define the size of the block
disp(['Encoding Block Size: ', num2str(blockSize)]);
ExpandedFrame = zeros(size(PadBarcode, 1) * blockSize, size(PadBarcode, 2) * blockSize, 3);

for channel = 1:3
    for i = 1:size(PadBarcode, 1)
        for j = 1:size(PadBarcode, 2)
            ExpandedFrame((i-1)*blockSize+1:i*blockSize, (j-1)*blockSize+1:j*blockSize, channel) = PadBarcode(i, j, channel);
        end
    end
end

% Save frame
imshow(ExpandedFrame, 'Border', 'loose');
drawnow;
pause(0.1); % Add a small pause to allow the graphics to update
F = getframe(gcf);
scaleFactor = 15; % Increase the scaleFactor to make the barcode cells larger and less blurry
ResizedFrame = imresize(F.cdata, scaleFactor);

ResizedFrame = imresize(ResizedFrame, 0.5);

% Ensure frame dimensions are even
if mod(size(ResizedFrame, 1), 2) ~= 0
    ResizedFrame = ResizedFrame(1:end-1,:,:);
end
if mod(size(ResizedFrame, 2), 2) ~= 0
    ResizedFrame = ResizedFrame(:,1:end-1,:);
end

[height, width, numChannels] = size(ResizedFrame);
disp(['Height: ', num2str(height)]);
disp(['Width: ', num2str(width)]);
disp(['Number of Channels: ', num2str(numChannels)]);

% Calculate new dimensions
new_height = ceil(height / 30) * 30;
new_width = ceil(width / 30) * 30;

% Resize the frame
ResizedFrame = imresize(ResizedFrame, [new_height new_width]);

% Check the new dimensions
[height, width, numChannels] = size(ResizedFrame);
disp(['New Height: ', num2str(height)]);
disp(['New Width: ', num2str(width)]);
disp(['Number of Channels: ', num2str(numChannels)]);

 % Save the resized frame as an image
    resized_filename = fullfile('frames', sprintf('resized_frame%d.png', frame));
    imwrite(ResizedFrame, resized_filename, 'png');


writeVideo(VideoUpload, ResizedFrame);

end

close(VideoUpload);
fclose('all');