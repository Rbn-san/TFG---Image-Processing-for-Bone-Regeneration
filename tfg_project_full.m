clear, clc, close all;

% ---------------------------- PARAMETERS ---------------------------------

LOW = 10;
MID = 15;
HIGH = 20;

% --- Sensitivity options:
% LOW : Better for maintain the shape of the elements
% MID : Compromise between high and low sensitivities
% HIGH : Better to differentiate tissues (recommended for images with
% osteoid tissue)

SENSITIVITY = LOW;

% ---------------------- Image Seleciton ----------------------------------

ImageName = uigetfile('*.jpg');
ImageRGB = imread(['images\selection\', ImageName]); % Change to your images folder

% -------------------- Color Deconvolution --------------------------------

    % INPUT PARAMETERS:
    ImgRGB                  = double(ImageRGB);
    StainingVectorID        = 1; % H&E
    DyeToBeRemovedID        = 0;
    doIcross                = 1;   

    % IMAGE RESHAPING:
    ImgR = ImgRGB(:,:,1);
    ImgG = ImgRGB(:,:,2);
    ImgB = ImgRGB(:,:,3);
    
    % CALL TO THE MAIN FUNCTION:
    [ImgR_back, ImgG_back, ImgB_back, Dye01_transmittance, Dye02_transmittance, Dye03_transmittance, LUTdye01, LUTdye02, LUTdye03, Q3x3Mat] = Colour_Deconvolution2(ImgR, ImgG, ImgB, StainingVectorID, DyeToBeRemovedID, doIcross);
    Img = uint8(Dye02_transmittance); % Eosin stain

% --- Dimensions
[M,N] = size(Img);

% ---------------------------- Markers ------------------------------------

ImgSelTitle = 'Please select the smallest trabeculae to be detected';
figure('Name', ImgSelTitle,'NumberTitle','off'); imshow(ImageRGB);
ROI = drawassisted;
ManualMarker = createMask(ROI);
roiArea = sum(sum(ManualMarker));
TRAB_SIZE = roiArea*0.6;

% --------------------------- TRABECULAE ----------------------------------

% --- Activity
seAct = strel("disk", 10);
ImgAct = imerode(imdilate(Img, seAct),seAct) - imdilate(imerode(Img, seAct),seAct);

% --- ACT_LVL and GRAY_LVL Calculation
grayMeanValue = round(sum(sum(ManualMarker.*double(Img)))/roiArea);
GRAY_LVL = grayMeanValue;

actMeanValue = round(sum(sum(ManualMarker.*double(ImgAct)))/roiArea);
ACT_LVL = actMeanValue + 5;

% --- Pixel Selection
PixSelTrab = Img > 5 & Img < GRAY_LVL + 100 & ImgAct < ACT_LVL;

% --- Labeling (con. 4)
ImgLabelTrab = im2uint8(bwareafilt(PixSelTrab, [TRAB_SIZE M*N], 4));

% --- Closing
seCl = strel("disk", SENSITIVITY);
ImgCl = imerode(imdilate(Img, seCl), seCl);

% --- Closing by Reconstruction
ImgRec =  imcomplement(imreconstruct(ImgLabelTrab, imcomplement(ImgCl)));

% --- Pixel Selection 2
PixSel2 = im2uint8(abs(int8(Img) - int8(ImgRec)) < 60 & ImgRec < GRAY_LVL + 50);

% --- Opening
seOp = strel("disk", 5);
ImgOp = imdilate(imerode(PixSel2, seOp), seOp);

% --- Hole Filling
MaskTrab = imfill(ImgOp, 'holes');

% --- Mask Application (RGB)
LogicMaskTrab = logical(MaskTrab);
ImgTrab = bsxfun(@times, ImageRGB, cast(LogicMaskTrab, 'like', ImageRGB));

% -------------------------- OSTEOID TISSUE -------------------------------

% --- White Balls detection -----------------------------------------------

% --- Pixel Selection
PixSel = Img > 250;

% --- Erotion
seEr = strel("square", 5);
ImgEroClearElem = imerode(PixSel, seEr);
 
% --- Labeling (con. 4)
ImgLabelWB = im2uint8(bwareafilt(ImgEroClearElem, [100 30000], 4));

% --- Trabecular elements deletion
ImgWB = ImgLabelWB - MaskTrab;

% --- Dilation + Erotion
seDil = strel("disk", 70); 
seEro = strel("disk", 50); 
MaskWB = imerode(imdilate(ImgWB, seDil), seEro);

% --- Osteoid tissue detection --------------------------------------------

% --- Pixel Selection
PixSel = im2uint8(Img < GRAY_LVL + 100 & ImgAct > ACT_LVL) - (MaskWB + MaskTrab);

% --- Erotion
seEr = strel("square", 5);
ImgEroOT = imerode(PixSel, seEr);

% --- Labeling (con. 4)
ImgLabelOT = im2uint8(bwareafilt(logical(ImgEroOT), [30000 M*N], 4));

% --- Closing
seCl = strel("disk", 20);
ImgClo = imerode(imdilate(ImgLabelOT, seCl),seCl);

% --- Hole Filling
MaskOT = imfill(ImgClo, 'holes') - MaskTrab;

% --- Mask Application (RGB)
LogicMaskOT = logical(MaskOT);
ImgOT = bsxfun(@times, ImageRGB, cast(LogicMaskOT, 'like', ImageRGB));

% ---------------------------- EMPTY SPACES -------------------------------

% --- Labeling (con. 4) Find white empty spaces
ImgLabelES = im2uint8(bwareafilt(ImgEroClearElem, [30000 M*N], 4));

% --- Hole Filling
MaskES = imfill(ImgLabelES, 'holes');

% ---------------------------- BONE MARROW --------------------------------

% --- Deletion of all elements that are not bone marrow
BM = imcomplement(MaskES + MaskTrab + MaskOT); 

% --- Noise removal (1. Opening - 2. Closing)
seNr = strel("disk", 8);
MaskBM = imerode(imdilate(imdilate(imerode(BM, seNr), seNr), seNr), seNr);

% --- Mask Application (RGB)
LogicMaskBM = logical(MaskBM);
ImgBM = bsxfun(@times, ImageRGB, cast(LogicMaskBM, 'like', ImageRGB));

% -------------------------- EVERYTHING ELSE ------------------------------

% --- Mask Application (RGB)
LogicMaskEE = logical(MaskES);
ImgEE = bsxfun(@times, ImageRGB, cast(LogicMaskEE, 'like', ImageRGB));

% --- Is every pixel classified? (Optional Check)
totalMask = MaskTrab + MaskOT + MaskBM + MaskES;

% -------------------------------------------------------------------------
% ------------------------- IMAGE CORRECTION ------------------------------

% Show results and ask for correction
% Order to follow when asking: 
% 1. Trabeculae 
% 2. Osteoid Tissue 
% 3. Bone Marrow
% 4. Everything Else

ImgTrabTitle = 'Trabeculae. Select: 1 Osteoid Tissue, 2 Bone Marrow, 3 Others';
figure('Name', ImgTrabTitle,'NumberTitle','off'); imshow(ImgTrab);
TrabROI_OT = drawassisted;
TrabROI_BM = drawassisted;
TrabROI_EE = drawassisted;

ImgOtTitle = 'Osteoid Tissue. Select: 1 Trabeculae, 2 Bone Marrow, 3 Others';
figure('Name', ImgOtTitle,'NumberTitle','off'); imshow(ImgOT);
OtROI_Trab = drawassisted;
OtROI_BM = drawassisted;
OtROI_EE = drawassisted;

ImgBmTitle = 'Bone Marrow. Select: 1 Trabeculae, 2 Osteoid Tissue, 3 Others';
figure('Name', ImgBmTitle,'NumberTitle','off'); imshow(ImgBM);
BmROI_Trab = drawassisted;
BmROI_OT = drawassisted;
BmROI_EE = drawassisted;

ImgEeTitle = 'Others. Select: 1 Trabeculae, 2 Osteoid Tissue, 3 Bone Marrow';
figure('Name', ImgEeTitle,'NumberTitle','off'); imshow(ImgEE);
EeROI_Trab = drawassisted;
EeROI_OT = drawassisted;
EeROI_BM = drawassisted;

% Check if there's something to correct
if ~(isempty(TrabROI_OT.Position) && isempty(TrabROI_BM.Position) && isempty(TrabROI_EE.Position) && ...
     isempty(OtROI_Trab.Position) && isempty(OtROI_BM.Position) && isempty(OtROI_EE.Position) && ...
     isempty(BmROI_Trab.Position) && isempty(BmROI_OT.Position) && isempty(BmROI_EE.Position) && ...
     isempty(EeROI_Trab.Position) && isempty(EeROI_OT.Position) && isempty(EeROI_BM.Position))
    
    Mask_TrabROI_OT = createMask(TrabROI_OT);
    Mask_TrabROI_BM = createMask(TrabROI_BM);
    Mask_TrabROI_EE = createMask(TrabROI_EE);

    Mask_OtROI_Trab = createMask(OtROI_Trab);
    Mask_OtROI_BM = createMask(OtROI_BM);
    Mask_OtROI_EE = createMask(OtROI_EE);

    Mask_BmROI_Trab = createMask(BmROI_Trab);
    Mask_BmROI_OT = createMask(BmROI_OT);
    Mask_BmROI_EE = createMask(BmROI_EE);

    Mask_EeROI_Trab = createMask(EeROI_Trab);
    Mask_EeROI_OT = createMask(EeROI_OT);
    Mask_EeROI_BM = createMask(EeROI_BM);

    % --- First of all, calculate the ROI masks to be fixed

    % Trabeculae
    Barrier = Mask_TrabROI_OT | Mask_TrabROI_BM | Mask_TrabROI_EE;
    AdditionMarkerTrab = Mask_OtROI_Trab | Mask_BmROI_Trab | Mask_EeROI_Trab;

    % Osteoid Tissue
    ExclusionMaskOT = im2uint8((Mask_OtROI_Trab | Mask_OtROI_BM | Mask_OtROI_EE) & LogicMaskOT);
    AdditionMaskOT = im2uint8((Mask_TrabROI_OT & LogicMaskTrab) | (Mask_BmROI_OT & LogicMaskBM) | (Mask_EeROI_OT & LogicMaskEE));

    % Bone Marrow
    ExclusionMaskBM = im2uint8((Mask_BmROI_Trab | Mask_BmROI_OT | Mask_BmROI_EE) & LogicMaskBM);
    AdditionMaskBM = im2uint8((Mask_TrabROI_BM & LogicMaskTrab) | (Mask_OtROI_BM & LogicMaskOT) | (Mask_EeROI_BM & LogicMaskEE));

    % Everything Else
    ExclusionMaskEE = im2uint8((Mask_EeROI_Trab | Mask_EeROI_OT | Mask_EeROI_BM) & LogicMaskEE);
    AdditionMaskEE = im2uint8((Mask_TrabROI_EE & LogicMaskTrab) | (Mask_OtROI_EE & LogicMaskOT) | (Mask_BmROI_EE & LogicMaskBM));

    % --- Trabeculae correction

    % Barrier
    ImgWithBarrier = ImgCl + im2uint8(Barrier * 255);
    ImgRecB = imcomplement(imreconstruct(ImgLabelTrab, imcomplement(ImgWithBarrier)));
    PixSelB = im2uint8(abs(int8(Img) - int8(ImgRecB)) < 60 & ImgRecB < GRAY_LVL + 50);
    
    % Addition
    TrabRoiArea = sum(sum(AdditionMarkerTrab));
    if TrabRoiArea == 0 
        grayMeanValue2 = GRAY_LVL;
    else
        grayMeanValue2 = round(sum(sum(AdditionMarkerTrab.*double(Img)))/TrabRoiArea);
    end
    seEr = strel("square", 5);                        
    ImgEroTCM = imerode(AdditionMarkerTrab, seEr);  
    ImgRecA =  imcomplement(imreconstruct(im2uint8(ImgEroTCM), imcomplement(ImgWithBarrier)));
    PixSelA = im2uint8(abs(int8(Img) - int8(ImgRecA)) < 60 & ImgRecA < grayMeanValue2 + 50);

    % Correction
    PixSelC = PixSelB + im2uint8(logical(PixSelA) & AdditionMarkerTrab);      
    ImgOp = imdilate(imerode(PixSelC, seOp), seOp);
    MaskTrab = imfill(ImgOp, 'holes');
    LogicMaskTrab = logical(MaskTrab);
    ImgTrab = bsxfun(@times, ImageRGB, cast(LogicMaskTrab, 'like', ImageRGB));

    % --- Osteoid Tissue correction
    MaskOT =  MaskOT + AdditionMaskOT - (MaskTrab + ExclusionMaskOT);
    LogicMaskOT = logical(MaskOT);
    ImgOT = bsxfun(@times, ImageRGB, cast(LogicMaskOT, 'like', ImageRGB));

    % --- Bone Marrow correction
    MaskBM = MaskBM + AdditionMaskBM - (MaskTrab + MaskOT + ExclusionMaskBM);
    LogicMaskBM = logical(MaskBM);
    ImgBM = bsxfun(@times, ImageRGB, cast(LogicMaskBM, 'like', ImageRGB));

    % --- "Everything Else" correction
    MaskEE = MaskES + AdditionMaskEE - (MaskTrab + MaskOT + MaskBM + ExclusionMaskEE);
    MaskEE = MaskEE + imcomplement(MaskTrab + MaskOT + MaskBM + MaskEE);
    LogicMaskEE = logical(MaskEE);
    ImgEE = bsxfun(@times, ImageRGB, cast(LogicMaskEE, 'like', ImageRGB));

end

% -------------------------------------------------------------------------
% --------------------------- FINAL DISPLAYS ------------------------------

figure; imshow(ImageRGB);

% --- Trabeculae
figure; imshow(ImgTrab);

% --- Osteoid Tissue
figure; imshow(ImgOT);

% --- Bone Marrow
figure; imshow(ImgBM);

% --- Everything Else
figure; imshow(ImgEE);

% -------------------------------------------------------------------------
% ------------------------- BONUS CALCULATIONS ----------------------------

% --- Osteocytes Counter
seOC = strel("square", 15);
TrabMaskOC = imerode(LogicMaskTrab, seOC);
TrabOC = uint8(TrabMaskOC).*Img;
PixSelOC = TrabOC > GRAY_LVL + 20 & ImgAct > 40;
Osteocytes = bwareafilt(PixSelOC,[1 100]);
[~,numOC] = bwlabel(Osteocytes);
figure; imshow(Osteocytes);

% --- Osteoblasts Area
% seOB1 = strel("disk", 23);
% seOB2 = strel("disk", 17);
% TrabMaskOB = imdilate(LogicMaskTrab, seOB1) - imerode(LogicMaskTrab, seOB2);
% TrabOB = uint8(TrabMaskOB).*Img;
% figure; imshow(TrabOB);

% --- OT and BM areas (in pixels)
pixNumOT = sum(sum(LogicMaskOT));
pixNumBM = sum(sum(LogicMaskBM));
