tic;
% 1. Read the video
VideoUpload = VideoReader('test1.mp4');

% Preallocate an array to hold the barcode frames
numFrames = VideoUpload.NumFrames;
FrameNumber= cell(numFrames, 1);
BarcodeFrames= cell(numFrames, 1);
referenceColorCellsArray = cell(1, numFrames);  
binaryReferenceArray = cell(1, numFrames);



% Set structuring element for morphological operations
se = strel('square', 3); % Modify size and shape according to your needs

previousLargestArea = 0;  % initialize to zero for the first frame
tolerance = 0.1;  % allow a 20% difference


% 2. Process each frame
for k = 1 : numFrames
    % Read the frame
    frame = read(VideoUpload, k);
    
    % Convert to grayscale to make it easier to threshold
    grayFrame = rgb2gray(frame);

    % Adjust the threshold level as needed
    thresholdLevel = 0.3;  % Adjust this value accordingly
    binaryFrame = imbinarize(grayFrame, thresholdLevel);
    binaryFrame = ~binaryFrame;
    
    binaryFrame = imerode(binaryFrame, strel('square',5));
   
    % Label the regions in the frame
    label = bwlabel(binaryFrame);

    % Find largest area to find position of barcode
    stats = regionprops(label,'Area');
    areaArray = [stats.Area];
    [sortedArea, sortedIndices] = sort(areaArray, 'descend');
    
    % If this is the first frame, use the largest area
    if k == 1 
        index = sortedIndices(1);
        previousLargestArea = sortedArea(1);
    else
        % Otherwise, find the largest area that is not significantly larger than the previous one
        for i = 1:length(sortedArea)
            if abs(sortedArea(i) - previousLargestArea) / previousLargestArea < tolerance
                index = sortedIndices(i);
                previousLargestArea = sortedArea(i);
                break;
            end
        end
    end
    
    BorderArea = (label==index);
    FilledImage = imfill(BorderArea,'holes');

    % Final barcode
    b = binaryFrame.*FilledImage;
    % Detect corners
    temp = FilledImage;
    se = strel('square',5);
    temp = imopen(temp,se);
    temp = imclose(temp,se);
    vertex = corner(temp,4);

    % Determine mid-point (origin) coordinates for the 4 quadrants
    mid_width = (max(vertex(:,1))+ min(vertex(:,1))) / 2;
    mid_height = (max(vertex(:,2))+ min(vertex(:,2))) / 2;

    % Identify the 4 quadrants
    LHS = vertex(vertex(:,1)<mid_width,:);
    RHS = vertex(vertex(:,1)>mid_width,:);
    Q1 = RHS(RHS(:,2)<mid_height,:); % Top Right
    Q2 = LHS(LHS(:,2)<mid_height,:); % Top Left
    Q3 = LHS(LHS(:,2)>mid_height,:); % Bottom Left
    Q4 = RHS(RHS(:,2)>mid_height,:); % Bottom Right

    % Projective transformation
    U = [Q1; Q2; Q3; Q4];
    T = fitgeotrans(U,[350 1;1 1;1 350;350 350],'projective');
    I_transformed = imwarp(frame,T,'OutputView',imref2d([350 350]));
    % Save the corrected barcode frame for later decoding
    BarcodeFrames{k} = I_transformed;
   
toc;

    % Get the current barcode frame
    barcodeFrame = BarcodeFrames{k};
    % Convert to grayscale
    grayFrame = rgb2gray(barcodeFrame);
    % Threshold the image to identify the white areas
    % Adjust the threshold level as needed
    thresholdLevel = 0.7;  % Adjust this value accordingly
    binaryFrame = imbinarize(grayFrame, thresholdLevel);
    % Find the regions
    stats = regionprops('table', binaryFrame, 'Area', 'BoundingBox');
    
    % Now, 'stats' contains information about each white region in the frame.
    % You can find the two largest regions, which should correspond to the barcode sections.
    [~, sortedIndices] = sort(stats.Area, 'descend');
    
    % Extract the bounding boxes of the two largest regions
    boundingBox1 = stats.BoundingBox(sortedIndices(1), :);
    boundingBox2 = stats.BoundingBox(sortedIndices(2), :);
    
    % Extract these areas
    section1 = imcrop(barcodeFrame, boundingBox1);
    section2 = imcrop(barcodeFrame, boundingBox2);
    
     % Further process to remove white border
    section1 = remove_white_border(section1, thresholdLevel);
    BarcodeFrames{k} =section1;
    
  toc;
headerRow = section2;

% Define the block size
blockSize = [10 10];  % 10x10 pixel blocks to match your cells

% Create a function to calculate the average color of a block
averageColorFcn = @(block_struct) mean(mean(block_struct.data, 1), 2);

% Apply the function to each block
% 'averageColorFcn' will be applied to each color channel separately
averageColorHeaderRow = blockproc(headerRow, blockSize, averageColorFcn);
% Remove the first and last rows
averageColorHeaderRow = averageColorHeaderRow(2:end-1, :, :);
% Remove the first and last columns
averageColorHeaderRow = averageColorHeaderRow(:, 2:end-1, :);
% Extract the reference color cells from the header image
referenceColorCells = averageColorHeaderRow(:, 22:29, :);
% Initialize a cell array to store the binary representations
binaryReference = cell(1, size(referenceColorCells, 2));

% Assign each reference color a binary value
for i = 1:size(referenceColorCells, 2)
    switch i
        case 1
            binaryReference{i} = [0 0 0];
        case 2
            binaryReference{i} = [0 0 1];
        case 3
            binaryReference{i} = [0 1 0];
        case 4
            binaryReference{i} = [0 1 1];
        case 5
            binaryReference{i} = [1 0 0];
        case 6
            binaryReference{i} = [1 0 1];
        case 7
            binaryReference{i} = [1 1 0];
        case 8
            binaryReference{i} = [1 1 1];
    end
end

% Define a function to compute the Euclidean distance between two colors
colorDistance = @(color1, color2) sqrt(sum((color1 - color2).^2));

% Initialize an array to store the binary representation of the header
headerBinary = zeros(size(averageColorHeaderRow, 2), 3);
% Iterate over each cell in the header
for i = 1:size(averageColorHeaderRow, 2)
    % Get the color of the current cell
    cellColor = averageColorHeaderRow(:, i, :);
    
    % Initialize a variable to store the minimum distance found
    minDistance = inf;
    % Initialize a variable to store the index of the closest color
    closestColorIndex = 0;

    % Iterate over each reference color
    for j = 1:size(referenceColorCells, 2)
        % Get the color of the current reference cell
        referenceColor = referenceColorCells(:, j, :);

        % Compute the distance between the cell color and the reference color
        distance = colorDistance(cellColor, referenceColor);

        % If this distance is smaller than the current minimum distance
        if distance < minDistance
            % Update the minimum distance
            minDistance = distance;
            % Update the index of the closest color
            closestColorIndex = j;
        end
    end

    % Assign the binary value of the closest color to the current cell
    headerBinary(i, :) = binaryReference{closestColorIndex};
end

% Transpose each row and then linearize
headerBinaryVec = reshape(headerBinary.', 1, []);

% Now you can extract the bits
errorBits = headerBinaryVec(1:3);
frameNumberBits = headerBinaryVec(4:12);


if isequal(errorBits, [1 0 0])
    disp('Error bits are [1 0 0]')
elseif isequal(errorBits, [0 1 0])
    disp('Error bits are [0 1 0]')
else
    disp('Frame should be eliminated')
end
frameNumberStr = num2str(frameNumberBits(:)');
frameNumber = bin2dec(frameNumberStr);

TotalFrameNum = headerBinaryVec(13:21);
TotalFrameNumStr = num2str(TotalFrameNum(:)');
TotalFrame = bin2dec(TotalFrameNumStr);
% Check if frameNumber is higher than TotalFrame
    if frameNumber > TotalFrame
        % If it is, skip this iteration and continue with the next frame
        continue
    end
FrameNumber{k} = frameNumber;
referenceColorCellsArray{k} = referenceColorCells;
binaryReferenceArray{k} = binaryReference;
toc;
end
toc;
% Create new arrays to hold the unique frames and their corresponding data
uniqueFrames = {};
uniqueReferenceColorArray = {};
uniqueBinaryReferenceArray = {};
uniqueFrameNumber = {};

% Initialize variables
lastSequenceStart = FrameNumber{1};
startOfSequence = 1;

% This variable will track the starting frame number of the last sequence we kept
lastKeptSequenceStart = -1;
seenFrameNumbers = [];  % list of all frame numbers that have been encountered

% Loop over all frames
for k = 2:length(FrameNumber)
    % If the frame number has changed, we've reached the end of a sequence
    if FrameNumber{k} ~= lastSequenceStart
        % Only keep this sequence if its starting frame number is different from the last one we kept
        % and it has not been seen before
        if lastSequenceStart ~= lastKeptSequenceStart && ~ismember(lastSequenceStart, seenFrameNumbers)
            % Select the middle frame of the sequence
            middleOfSequence = startOfSequence + floor((k - startOfSequence) / 2);

            % Save the middle frame and its corresponding data
            uniqueFrames{end+1} = BarcodeFrames{middleOfSequence};
            uniqueReferenceColorArray{end+1} = referenceColorCellsArray{middleOfSequence};
            uniqueBinaryReferenceArray{end+1} = binaryReferenceArray{middleOfSequence};
            uniqueFrameNumber{end+1} = FrameNumber{middleOfSequence};

            % Update the last kept sequence start
            lastKeptSequenceStart = lastSequenceStart;

            % Add the frame number to the seen list
            seenFrameNumbers = [seenFrameNumbers, lastSequenceStart];
        end
        
        % Update variables for the next sequence
        lastSequenceStart = FrameNumber{k};
        startOfSequence = k;
    end
end

% Handle the last sequence, which isn't handled by the loop
if lastSequenceStart ~= lastKeptSequenceStart && ~ismember(lastSequenceStart, seenFrameNumbers)
    middleOfSequence = startOfSequence + floor((length(FrameNumber) - startOfSequence) / 2);
    uniqueFrames{end+1} = BarcodeFrames{middleOfSequence};
    uniqueReferenceColorArray{end+1} = referenceColorCellsArray{middleOfSequence};
    uniqueBinaryReferenceArray{end+1} = binaryReferenceArray{middleOfSequence};
    uniqueFrameNumber{end+1} = FrameNumber{middleOfSequence};
end

% Sort FrameNumber in ascending order
[sortedFrameNumber, sortOrder] = sort([uniqueFrameNumber{:}]);

% Reorder the other arrays to match
uniqueFrames = uniqueFrames(sortOrder);
uniqueReferenceColorArray = uniqueReferenceColorArray(sortOrder);
uniqueBinaryReferenceArray = uniqueBinaryReferenceArray(sortOrder);

% Convert sortedFrameNumber back to cell array
uniqueFrameNumber = num2cell(sortedFrameNumber);

toc;
% Define a function to compute the Euclidean distance between two colors
colorDistance = @(color1, color2) sqrt(sum((color1 - color2).^2));

% 3. Decode the barcode data 
% Iterate over each unique frame
FrameBinary = [];
for k = 1:length(uniqueFrames)
    Framedata = uniqueFrames{k};
    FrameCells = uniqueReferenceColorArray{k};
    FrameBinaryReference = uniqueBinaryReferenceArray{k};
    % Get the size of the frame data
    [frameHeight, frameWidth, ~] = size(Framedata);

    % Calculate the width and height of a cell
    cellWidth = frameWidth / 30;
    cellHeight = frameHeight / 30;
H = ones(3,3) / 9;  % 3x3 mean filter
    % Iterate over each cell
    for i = 1:30
        for j = 1:30
            % Calculate the coordinates of the cell
            cellX = round((j - 1) * cellWidth) + 1;
            cellY = round((i - 1) * cellHeight) + 1;
            cellWidthActual = round(j * cellWidth) - cellX + 1;
            cellHeightActual = round(i * cellHeight) - cellY + 1;

            % Extract the cell from the frame data
             cellData = Framedata(max(cellY+1,1):min(cellY+cellHeightActual-2,frameHeight), max(cellX+1,1):min(cellX+cellWidthActual-2,frameWidth), :);
             % Apply the mean filter to the cellData
        for l = 1:3
            cellData(:,:,l) = filter2(H, double(cellData(:,:,l)));
        end
            % Calculate the average color of the cell
            averageColor = mean(mean(cellData, 1), 2);

            % Initialize a variable to store the minimum distance found
            minDistance = inf;
            % Initialize a variable to store the index of the closest color
            closestColorIndex = 0;

            % Iterate over each reference color
            for n = 1:size(FrameCells, 2)
                % Get the color of the current reference cell
                referenceColor = FrameCells(:, n, :);

                % Compute the distance between the cell color and the reference color
                distance = colorDistance(averageColor, referenceColor);

                % If this distance is smaller than the current minimum distance
                if distance < minDistance
                    % Update the minimum distance
                    minDistance = distance;
                    % Update the index of the closest color
                    closestColorIndex = n;
                end
            end

            % Assign the binary value of the closest color to the current cell
            binaryValue = FrameBinaryReference{closestColorIndex};
            FrameBinary = [FrameBinary, binaryValue];  
            
        end
    end
end
toc;
   % 4. Convert the binary to ASCII or words
% Compute how many bits we are short of a multiple of 14
paddingLength = mod(length(FrameBinary), 14);
if paddingLength ~= 0
    paddingLength = 14 - paddingLength;  % The actual number of bits needed
end

% Add the padding to the end of FrameBinary
FrameBinary = [FrameBinary, zeros(1, paddingLength)];

% Reshape FrameBinary into a matrix where each row represents a 14-bit character
binaryMatrix = reshape(FrameBinary, 14, [])';

% Initialize a string to store the decoded message
decodedMessage = '';

% Remove the first 7 columns (unwanted bits)
binaryMatrix = binaryMatrix(:, 8:end);

% Loop over each row in the binaryMatrix
for i = 1:size(binaryMatrix, 1)
    % Get the binary code for the current character
    binaryCode = binaryMatrix(i, :);
    
    % Convert the binary code into a decimal number
    decimalNumber = bin2dec(sprintf('%1d', binaryCode));
    
    % Convert the decimal number into an ASCII character
    asciiChar = char(decimalNumber);
    
    % Append the ASCII character to the decoded message
    decodedMessage = [decodedMessage, asciiChar];
end

fileID = fopen('decodedMessage.txt', 'w');  % 'w' stands for write mode
if fileID == -1
    error('Cannot open file for writing.');
end
fprintf(fileID, '%s', decodedMessage);
fclose(fileID);

% Close the video file
clear VideoUpload;
toc;
% Function to remove white border from an image
function img_no_border = remove_white_border(img, thresholdLevel)
    binaryImg = imbinarize(rgb2gray(img), thresholdLevel);
    binaryImg = ~binaryImg;  % Invert binary image to make white pixels as background
    stats = regionprops('table', binaryImg, 'Area', 'BoundingBox');
    [~, index] = max(stats.Area);
    boundingBox = stats.BoundingBox(index, :);
    img_no_border = imcrop(img, boundingBox);
end