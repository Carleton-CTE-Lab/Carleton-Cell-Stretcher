%% ========================================================================
% Script Name: PI_Cell_Quantification.m
%  Author: Gia Kang
%  Institution: CTE Lab, Carleton University, Canada
%  Date: Oct 2025
%  Version: 1.0
%
% Purpose:
%   This script quantifies membrane-compromised (PI-positive) cells from 
%   fluorescence microscopy images of nuclei (Hoechst/DAPI channel) and 
%   propidium iodide (PI) staining. It compares pixel intensity in the PI 
%   channel to background intensity to identify damaged cells.
%
% Description:
%   The script processes paired fluorescence images stored as .tif files 
%   (named according to the convention: "FileName_1.tif" for Hoechst and 
%   "FileName_2.tif" for PI). For each image pair, it:
%       1. Loads and normalizes the images to 16-bit depth.
%       2. Generates binary masks from the nuclei channel.
%       3. Removes large or small artifacts based on area thresholds.
%       4. Measures mean PI intensity for each nucleus.
%       5. Classifies a nucleus as PI-positive if intensity exceeds a
%          user-defined threshold (bkINT * threshold).
%       6. Saves summary results to an Excel file.
%
% Input:
%   - Folder containing .tif image pairs of nuclei and PI stains.
%
% Output:
%   - Excel file with results:
%       Column 1: Total nuclei count
%       Column 2: PI-positive nuclei count
%       Column 3: % PI-positive cells
%
% Parameters to Adjust:
%   - threshold: sets sensitivity for PI positivity (default = 1.25)
%   - depth: image bit depth (default = 65536 for 16-bit)
%   - Area thresholds in 'bwareaopen' may be adjusted based on cell size.
%
% Required Software:
%   - Image Processing Toolbox
%
% ========================================================================


%% INITIALIZATION
clear all; close all; % clear workspace and close all open figures

% Prompt user to select any .tif file from the dataset folder
[fname, fpath] = uigetfile('*.tif');

% Retrieve all .tif files from that folder
a = dir([fpath '*.tif']);
nfile = max(size(a)); % number of image files detected

% Sort file list by numerical order (ensures correct pairing, e.g. 1, 2, 3...)
[~, reindex] = sort(str2double(regexp({a.name}, '\d+', 'match', 'once')));
a = a(reindex);

names = {}; % store sample names for output labeling


%% USER-DEFINED PARAMETERS
depth = 65536;                 % image bit depth (16-bit = 65536 gray levels)
threshold = 1.5;              % intensity threshold multiplier for PI-positivity
saveFileName = 'Sample_Results'; % name of Excel results file
results = zeros(nfile/2, 3);   % preallocate results array


%% FILE NAMING CONVENTION
% Ensure each pair of files follows this naming format:
% "SampleName_1.tif" → Hoechst channel
% "SampleName_2.tif" → Propidium Iodide (PI) channel


%% MAIN ANALYSIS LOOP
for n = 1:nfile/2

    % Get file names for Hoechst and PI channels
    fname_hc = a(n*2-1).name;
    fname_pi = a(n*2).name;

    % Extract sample name (remove trailing "_1.tif" or "_2.tif")
    names{end+1} = fname_hc(1:end-6);

    % Create figure for visualization of intermediate steps
    h = figure;

    %% Load and display Hoechst (DAPI) image
    if str2double(fname_hc(end-4)) == 1
        S1 = imread([fpath fname_hc]);
        I1 = mat2gray(S1, [0 depth]);  % normalize image intensity to [0,1]
        I1_rgb = ind2rgb(imadjust(S1),abyss(depth));    %colorized image for figure
        subplot(3,2,1);
        imshow(I1_rgb);
        title('Nuclei (DAPI) stain');
    end

    %% Load and display Propidium Iodide (PI) image
    if str2double(fname_pi(end-4)) == 2
        S2 = imread([fpath fname_pi]);
        I2 = mat2gray(S2, [0 depth]);  % normalize image intensity to [0,1]
        I2_rgb = ind2rgb(imadjust(S2),hot(depth*1.5));   %colorized image for figure
        subplot(3,2,2);
        imshow(I2_rgb);
        title('PI stain');
    end

    %% Create binary mask of nuclei (from DAPI channel)
    MASK = imbinarize(I1); % threshold-based binarization
    bkMASK = imcomplement(MASK); % invert mask to get background region

    %% Calculate average background intensity in PI image
    bkINT = mean2(I2(find(bkMASK))); % mean intensity of background pixels

    %% Remove unwanted artifacts from the nuclei mask
    bigBlobs = bwareaopen(MASK, 10000); % remove very large clumps (e.g., cell clusters)
    smallBlobs = MASK - bigBlobs;       % isolate small/medium nuclei
    smallBlobs = bwareaopen(smallBlobs, 100); % remove small debris
    areaSmall = bwarea(smallBlobs);     % total area of valid nuclei

    % Visualize cleaned nuclei mask
    subplot(3,2,3);
    imshow(smallBlobs);
    title('Cleaned nuclei mask');

    % Apply nuclei mask to PI channel to isolate nuclear regions
    subplot(3,2,4);
    imshow(I2_rgb .* smallBlobs); %colorized image for figure
    title('PI channel (masked by nuclei)');

    %% Measure PI intensity per nucleus
    stats = regionprops(smallBlobs, I2, 'MeanIntensity');
    results(n,1) = length(stats); % total number of detected nuclei

    % Extract mean intensity of PI signal for each nucleus
    for m = 1:length(stats)
        Intensity(m) = stats(m).MeanIntensity;
    end

    %% Generate histogram and visualize PI intensity distribution
    subplot(3,2,[5,6]);
    hold on;
    histogram(Intensity);
    plot(bkINT * threshold, 0, 'r*'); % red star indicates threshold cutoff
    title('PI Intensity Distribution');
    hold off;

    %% Identify PI-positive nuclei
    positives = find(Intensity > bkINT * threshold);

    % Store numerical results
    results(n,2) = length(positives);                      % number of PI-positive nuclei
    results(n,3) = (length(positives) / length(stats)) * 100; % percentage of PI-positive nuclei

    saveFigFolder = fullfile(fpath, 'AnalysisFigures');
    if ~exist(saveFigFolder, 'dir')
        mkdir(saveFigFolder);  % create folder if it doesn't exist
    end

    % After finishing all plotting for the iteration (before clearing)
    figFileName = fullfile(saveFigFolder, [names{n}, '_analysis.png']);  % filename per sample
    saveas(h, figFileName);  % saves figure as PNG
    close(h);

    % Optionally save as MATLAB figure (.fig) for later editing
    % figFileNameFig = fullfile(saveFigFolder, [names{n}, '_analysis.fig']);
    % savefig(h, figFileNameFig);

    %% Clean up temporary variables before next iteration
    clear positives S1 S2 I1 I2 MASK smallBlobs bigBlobs stats Intensity
    
end


%% OUTPUT RESULTS
%% Combine names and results for Excel output
% Convert numeric results to cell array
resultsCell = num2cell(results);

% Combine names as the first column
names = transpose(names);
outputCell = [names, resultsCell];

% Optionally add header row
header = {'SampleName', 'TotalNuclei', 'PI_Positive', '%PI_Positive'}; 
outputCell = [header; outputCell];

% Write combined data to Excel
excelFileName = fullfile(fpath, strcat(saveFileName, '.xlsx'));
writecell(outputCell, excelFileName);

fprintf("Analysis complete! Results saved as %s.xlsx\n", saveFileName);
