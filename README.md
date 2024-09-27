# 3D Colour Barcode for High-Capacity Data Transfer

# Project Overview

This project focuses on developing a robust 3D-coloured barcode system that enhances the data capacity of traditional barcodes by leveraging 3D imaging and colour encoding techniques. The barcode system was implemented using MATLAB and enables efficient encoding and decoding of large amounts of data, without relying on traditional network methods like Wi-Fi or Bluetooth. The project demonstrates the potential of 3D barcodes for offline data transfer in diverse applications such as advertising and education.

# Objective

  - High Data Capacity: Design a 3D barcode system capable of encoding and decoding large data volumes.
  - Offline Data Transfer: Develop a system that enables data transfer via a coloured video barcode displayed on a monitor and recorded with a mobile device.
  - Efficient Algorithms: Implement MATLAB-based algorithms for encoding, decoding, and error detection.
    
# Features

  - 3D Barcode Design: Utilizes a matrix of coloured cells to represent binary data, enhancing data capacity.
  - MATLAB Algorithms: Encoding and decoding algorithms implemented in MATLAB.
  - Error Detection: Basic error detection included for frame validation.
  - Offline Data Transfer: The barcode can be displayed and scanned without needing internet connectivity.

# Technologies Used

  - MATLAB: For coding the encoding and decoding algorithms.
  - 3D Imaging: The barcodes are displayed as a sequence of frames that simulate a 3D structure through time multiplexing.
  - Colour Encoding: A combination of 8 different colours (black, blue, cyan, green, white, red, yellow, magenta) is used to encode data efficiently.

# How to Use the Project

  1. Clone the Repository:
     ```bash
     git clone https://github.com/Puteraharis24/HighCapacity-3D-Barcode.git

  2. Encoding process:

     - Run the Encoding.m script in MATLAB to convert your text or data into a series of 3D barcode frames.
     - The barcode is displayed as a looping video of frames representing the encoded data.

  3. Decoding process:

     - Record the barcode video with a mobile device.
     - Run the Decoding.m script in MATLAB to decode the recorded video and extract the original data.

  4. Requirements:

     - MATLAB R2023 or later
     - A monitor to display the barcodes
     - A mobile device with a camera to record the barcode video
    
# Project workflow

  1. Encoding:

     - Data is converted into binary and divided into 30x30 cells.
     - Each cell is assigned a colour based on the binary data (3 bits per cell).
     - The encoded data is displayed a a series of frames.

  2. Recording:

     - The barcode frames are displayed on a monitor.
     - A mobile device records the video for later decoding.

  3. Decoding:

     - The video is processed to extract individual frames.
     - The MATLAB decoding algorithm converts the colours back to binary data.
     - The binary data is reassembled into the original data.

# Results

  - The system achieved aan 80-90% decoding success rate when recorded in an upright position. However, accuracy decreases as the recording angle increases due to parallax distortion.

# Future Work

  - Error Correction: Implement advanced error correction techniques for more robust data transfer.
  - Mobile App: Develop a mobile application for easier scanning and decoding.
  - Optimizations: Enhance decoding speed and accuracy under various lighting and angle conditions.

# License

This project is licensed under the MIT License - see the LICENSE file for details.

  
