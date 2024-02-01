%% Noisy Image analysis
%August 2022

close all; clearvars;

currpath = pwd;
Idx = regexp(pwd,'Images/','end');
main_path = currpath(1:Idx);
redoLocalPolys = 0;
%% 1) Load Experiment Data

fprintf('Loading Experiment Data\n');

redo_data_compile = 1; % 1 to rerun data compilations (compileDataExp1_3), 0 to use previous compilation  
if exist(fullfile(main_path, 'scripts/Analysis/Scoring.mat'), 'file')
    load(fullfile(main_path, 'scripts/Analysis/Scoring.mat'))
end

if ~exist('../../data/Exp1/AllData/AllExp1ControlData.mat', 'file') || ~exist('../../data/Exp1/AllData/AllExp1Data.mat', 'file') || redo_data_compile
    [AllExp1Data, AllControlData, AllCatchData] = compileDataExp1_3(currpath);
else
    load('../../data/Exp1/AllData/AllExp1ControlData.mat');
    load('../../data/Exp1/AllData/AllExp1Data.mat')
    load('../../data/Exp1/AllData/AllCatchData.mat');
end

fprintf('Loading Polygons\n');

load('images_withPoly.mat'); %exp polys (images_withPoly2 with new Crocodile Eye Poly)
load('CatchPolys.mat'); %catch polys
load('../LookupPromptsCatch.mat');
load('../LookupPromptsExp1.mat');


%% 2) Define Variables

fprintf('Defining Variables\n');

imageIDs = [1 2 4 5 7 8 11 12 13 14 15 16 17 24 25 26 27 32 39 40]; %IDs of images used in final study
NumIi = length(imageIDs);
ssIDs = unique(AllExp1Data.pertrial.pID);  %all participant IDs
NumSs = length(ssIDs);                     

sizeMark = 20;
pxpermm = 680/273.9; %680px = 273.9mm
mmperpx = 273.9/680;
redoPolys = 0; % 1 to redraw polygons do poiting accuracy, two-tone images
redoPolysCatch = 0;% 1 to redraw polygons do poiting accuracy, catch images

%colour maps
Rainbow = [0.265 0.555 0.465; 0.46 0.42 0.62; 0.87 0.485 0.6; 0.855 0.58 0.4];
Light = [0.455 0.745 0.655; 0.65 0.61 0.81; 1 0.675 0.795; 1 0.77 0.59];
MeanColours = {[0.8157 0.8157 0.8157], [0.6706  0.6706  0.6706], [0.5725 0.5725 0.5725], [0.1412 0.1412 0.1412]; [1 0.7490 0.247], [1 0.5529 0.4431], [1 0.4196 0.5765], [0.5098 0.2549 0.7529]}; 
spring2 = [linspace(0.4,0.6,66)',linspace(0.2,0.3,66)', linspace(0.7,0.8,66)'; ones(190,1),linspace(0,1,190)', linspace(1,0,190)' ]; %colourmap



x = jet;
jet2 = x(0.5*end:end,:);

BGsize = [1080 1920];
%% 3) Extract Per-Subject & Per-Image Data

if ~exist(fullfile(main_path, '/scripts/Analysis/Scoring.mat'), 'file')

    disp('Extracting Data');

    %define structs
    subs.ID = nan(NumSs,1);
    subs.ages = nan(NumSs,1);
    subs.AgeGroup = nan(NumSs,1);
    subs.Gender = nan(NumSs,1);

    subs.ImNum = nan(NumSs,NumIi);
    subs.TwoToneAcc = nan(NumSs,NumIi);
    subs.GrayScAcc = nan(NumSs,NumIi);
    subs.CanSeeNow = nan(NumSs,NumIi);

    subs.TwoToneCoord1 = nan(NumSs,NumIi,2);
    subs.TwoToneCoord2 = nan(NumSs,NumIi,2);
    subs.GrayScaleCoord1 = nan(NumSs,NumIi,2);
    subs.GrayScaleCoord2 = nan(NumSs,NumIi,2);

    for ss = 1:NumSs
        for ii = 1:NumIi
            ssIdx = ssIDs(ss);
            iiIdx = imageIDs(ii);
            idx = AllExp1Data.pertrial.pID == ssIdx & AllExp1Data.pertrial.ImNum == iiIdx;
            if sum(idx) == 1
                subs.ID(ss) = AllExp1Data.pertrial.pID(idx);
                subs.ages(ss) = AllExp1Data.pertrial.age(idx);
                subs.AgeGroup(ss) = AllExp1Data.pertrial.ageGroup(idx);

                subs.ImNum(ss,ii) = AllExp1Data.pertrial.ImNum(idx);
                subs.TwoToneReport(ss,ii) = AllExp1Data.pertrial.TwoToneReport(idx);
                subs.GreyScaleReport(ss,ii) = AllExp1Data.pertrial.GreyScaleReport(idx);
                subs.TwoToneAcc(ss,ii) = AllExp1Data.pertrial.TwoToneAcc(idx);
                subs.GrayScAcc(ss,ii) = AllExp1Data.pertrial.GrayScAcc(idx);
                subs.CanSeeNow(ss,ii) = AllExp1Data.pertrial.CanSeeNow(idx);

                subs.TwoToneCoord1(ss,ii,:) = AllExp1Data.pertrial.coord1(idx,:);
                subs.TwoToneCoord2(ss,ii,:) = AllExp1Data.pertrial.coord2(idx,:);
                subs.GrayScaleCoord1(ss,ii,:) = AllControlData.pertrial.coord1(idx,:);
                subs.GrayScaleCoord2(ss,ii,:) = AllControlData.pertrial.coord2(idx,:);


            elseif sum(idx) > 1 % if there are two of the same images in single particpant
                fprintf('subject %i',ss)
                fprintf('image %i',ii)
                break
            end
        end
    end
end


%Find Individual distances between TwoTone and GreyScale coordinates
%(pointing distance)
subs.TT2GSDist1 = mmperpx*hypot((subs.TwoToneCoord1(:,:,1) - subs.GrayScaleCoord1(:,:,1)), (subs.TwoToneCoord1(:,:,2) - subs.GrayScaleCoord1(:,:,2))); %distance in mm
subs.TT2GSDist2 = mmperpx*hypot((subs.TwoToneCoord2(:,:,1) - subs.GrayScaleCoord2(:,:,1)), (subs.TwoToneCoord2(:,:,2) - subs.GrayScaleCoord2(:,:,2))); %distance in mm

%% 3) Extract Per-Subject Catch Data

disp('Extracting Catch Data');

%define structs
subs.CaImNum = nan(NumSs,4);
subs.CaTwoToneAcc = nan(NumSs,4);
subs.CaGrayScAcc = nan(NumSs,4);
subs.CaCanSeeNow = nan(NumSs,4);

subs.CaTwoToneCoord1 = nan(NumSs,4,2);
subs.CaTwoToneCoord2 = nan(NumSs,4,2);
subs.CaGrayScaleCoord1 = nan(NumSs,4,2);
subs.CaGrayScaleCoord2 = nan(NumSs,4,2);

for ss = 1:NumSs
    for ii = 1:5
        newim = ii;
        if ii == 5
            newim = 4;
        end 

        ssIdx = ssIDs(ss);

        idx = AllCatchData.pertrial.pID == ssIdx & AllCatchData.pertrial.ImNum == ii;


        if sum(idx) == 1

                %exclude turtle images
    if contains(AllCatchData.pertrial.ImLabel(idx), 'Turtle')
        if ~cellfun(@isempty,AllCatchData.pertrial.GrayScReport{idx})
           if contains(AllCatchData.pertrial.GrayScReport{idx}, 'turtle') || contains(AllCatchData.pertrial.GrayScReport{idx}, 'tortoise')
        AllCatchData.pertrial.GrayScReport{idx}
        AllCatchData.pertrial.coord1(idx,1) = NaN; AllCatchData.pertrial.coord1(idx,2) = NaN;
        AllCatchData.pertrial.coord2(idx,1) = NaN; AllCatchData.pertrial.coord2(idx,2) = NaN;
        AllCatchData.pertrial.GrayScAcc(idx) = NaN; AllCatchData.pertrial.TwoToneAcc(idx) = NaN;
           end
        end
    end

            subs.CaImNum(ss,newim) = AllCatchData.pertrial.ImNum(idx);
            subs.CaTwoToneReport(ss,newim) = AllCatchData.pertrial.TwoToneReport(idx);
            subs.CaGreyScaleReport(ss,newim) = AllCatchData.pertrial.GrayScReport(idx);
            subs.CaTwoToneAcc(ss,newim) = AllCatchData.pertrial.TwoToneAcc(idx);
            subs.CaGrayScAcc(ss,newim) = AllCatchData.pertrial.GrayScAcc(idx);
            subs.CaCanSeeNow(ss,newim) = AllCatchData.pertrial.CanSeeNow(idx);

            subs.CaTwoToneCoord1(ss,newim,:) = AllCatchData.pertrial.coord1(idx,:);
            subs.CaTwoToneCoord2(ss,newim,:) = AllCatchData.pertrial.coord2(idx,:);
            subs.CaGrayScaleCoord1(ss,newim,:) = AllCatchData.pertrial.coord1(idx,:);
            subs.CaGrayScaleCoord2(ss,newim,:) = AllCatchData.pertrial.coord2(idx,:);


        elseif sum(idx) > 1 % if there are two of the same images in single particpant
            fprintf('subject %i',ss)
            fprintf('image %i',ii)
            break
        end
    end
end

subs.CaImNum = subs.CaImNum + 40;

%Find Individual distances between TwoTone and GreyScale coordinates
%(pointing distance)
subs.CaTT2GSDist1 = mmperpx*hypot((subs.CaTwoToneCoord1(:,:,1) - subs.CaGrayScaleCoord1(:,:,1)), (subs.CaTwoToneCoord1(:,:,2) - subs.CaGrayScaleCoord1(:,:,2))); %distance in mm
subs.CaTT2GSDist2 = mmperpx*hypot((subs.CaTwoToneCoord2(:,:,1) - subs.CaGrayScaleCoord2(:,:,1)), (subs.CaTwoToneCoord2(:,:,2) - subs.CaGrayScaleCoord2(:,:,2))); %distance in mm


%% Demographics

Ages.Groups = unique(subs.AgeGroup(~isnan(subs.AgeGroup)));
Ages.NumGroups = length(Ages.Groups);

for i = 1:Ages.NumGroups
    idx = subs.AgeGroup == i;
    Ages.mean(i) = mean(subs.ages(idx),'omitnan');
    Ages.minmax(i,1) = nanmin(subs.ages(idx)); %#ok<NANMIN> 
    Ages.minmax(i,2) = nanmax(subs.ages(idx)); %#ok<NANMAX> 
    Ages.SD(i) = std(subs.ages(idx),'omitnan');
    Ages.size(i) = sum(idx);
    Ages.GroupNames{i} = [num2str(floor(Ages.minmax(i,1))), '-', num2str(floor(Ages.minmax(i,2)))];
end

ColorRange = max(0.5*subs.ages); % for transferring colormap, - & + this number

%% Polygons
if redoPolys
    for ii = 1:NumIi %#ok<UNRCH> 
        iid = imageIDs(ii);
        background = NaN(1080,1920);

        GreyImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/Greyscale/Im_' num2str(iid), '.jpg']));
        GreyImage = imresize(GreyImage, [680 NaN]);
        GreyImageInsert = insertMatrix(background, GreyImage);

        imagesc(GreyImageInsert)
        colormap(gray)


        disp(LookupPrompt(iid).FirstPrompt)
        poly = drawpolygon('FaceAlpha',0,'Color','r');
        images.ROI1Poly{iid} = poly.Position;
        pause
        clearvars poly
        disp(LookupPrompt(iid).SecondPrompt)
        poly = drawpolygon('FaceAlpha',0,'Color','r');
        images.ROI2Poly{iid} = poly.Position;
        pause
        clearvars poly

    end

end

if redoLocalPolys
    for ii = 1:NumIi
        iid = imageIDs(ii);
        background = NaN(1080,1920);

        TwoToneImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);
        TwoToneInsert = insertMatrix(background, TwoToneImage);

        imagesc(TwoToneInsert)
        colormap(gray)

        disp(LookupPrompt(iid).FirstPrompt)
        poly = drawpolygon('FaceAlpha',0,'Color','r','Position',images.ROI1Poly{iid});
        pause
        images.ROI1LocalPoly{iid} = poly.Position

            bigpoly = scale(polyshape(images.ROI1Poly{iid}),sqrt(3),mean(images.ROI1Poly{iid}))
            bigpoly = drawpolygon('FaceAlpha',0,'Color','r','Position',bigpoly.Vertices);
            images.ROI1BigPoly{iid} = bigpoly.Position;
            pause

        clearvars poly
        disp(LookupPrompt(iid).SecondPrompt)
        poly = drawpolygon('FaceAlpha',0,'Color','r','Position',images.ROI2Poly{iid});
        pause
        images.ROI2LocalPoly{iid} = poly.Position

            bigpoly = scale(polyshape(images.ROI2Poly{iid}),sqrt(3),mean(images.ROI2Poly{iid}))
            bigpoly = drawpolygon('FaceAlpha',0,'Color','r','Position',bigpoly.Vertices);
            images.ROI2BigPoly{iid} = bigpoly.Position;
            pause
  
        clearvars poly

    end

end

if redoPolysCatch
    for ii = 1:length(CatchLookupPrompt) %#ok<UNRCH> 
        background = NaN(1080,1920);

        GreyImage = imread(fullfile(main_path, ['/stimuli/CatchTrialStimuli/Greyscale/Im_' num2str(ii), '.jpg']));
        GreyImage = imresize(GreyImage, [680 NaN]);
        GreyImageInsert = insertMatrix(background, GreyImage);


        imagesc(GreyImageInsert)
        colormap(gray)

        disp(CatchLookupPrompt(ii).FirstPrompt)
        poly = drawpolygon('FaceAlpha',0,'Color','r');
        Catch.ROI1Poly{ii} = poly.Position;
        pause
        clearvars poly
        disp(CatchLookupPrompt(ii).SecondPrompt)
        poly = drawpolygon('FaceAlpha',0,'Color','r');
        Catch.ROI2Poly{ii} = poly.Position;
        pause
        clearvars poly

    end
end


%plot polygons
figure
tiledlayout('flow', 'Padding', 'none')
for ii = 1:NumIi
    iid = imageIDs(ii);
    background = NaN(1080,1920);

    % load two tone (fig 17 weird border)
    TwoToneImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/TwoTone/Im_' num2str(iid), '.jpg']));
    TwoToneImage = imresize(TwoToneImage, [680 NaN]);
    TwoToneInsert = insertMatrix(background, TwoToneImage);

    % load greyscale
    GreyImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/Greyscale/Im_' num2str(iid), '.jpg']));
    GreyImage = imresize(GreyImage, [680 NaN]);
    GreyImageInsert = insertMatrix(background, GreyImage);


    ROI1 = polyshape(images.ROI1Poly{iid});
    ROI2 = polyshape(images.ROI2Poly{iid});
    
    if images.ROI1LocalPoly{iid} ~= images.ROI1Poly{iid}
        ROI1l = images.ROI1LocalPoly{iid};
    end

    if images.ROI2LocalPoly{iid} ~= images.ROI2Poly{iid}
        ROI2l = images.ROI2LocalPoly{iid};
    end
    set(gcf,'Color',[1,1,1])%%%set figure background to white

    nexttile

    imagesc(GreyImageInsert)
    colormap(gray)
    hold on

    plot(ROI1)
    hold on
    plot(ROI2)
    axis off

    nexttile

    imagesc(TwoToneInsert)
    colormap(gray)
    hold on

    plot(ROI1)
    hold on
    plot(ROI2)
    hold on
    plot(ROI1l)
    hold on
    plot(ROI2l)

    axis off
end

%% compute Polygon Accuracy
disp('Computing Catch Trial Polygon Accuracy');

%define variables
AllCatchData.pertrial.Coord1PolyAcc = nan(size(AllCatchData.pertrial.pID));
AllCatchData.pertrial.Coord2PolyAcc = nan(size(AllCatchData.pertrial.pID));
AllCatchData.pertrial.ROI1Area = nan(size(AllCatchData.pertrial.pID));
AllCatchData.pertrial.ROI2Area = nan(size(AllCatchData.pertrial.pID));

for tt = 1:length(AllCatchData.pertrial.pID)
    background = NaN(1080,1920);
    ii = AllCatchData.pertrial.ImNum(tt);

    ROI1 = polyshape(Catch.ROI1Poly{ii});
    ROI2 =  polyshape(Catch.ROI2Poly{ii});

    AllCatchData.pertrial.ROI1Area(tt) = area(ROI1);
    AllCatchData.pertrial.ROI2Area(tt) = area(ROI2);

    %Individual within bound?
    withinBound1 = inpolygon(AllCatchData.pertrial.coord1(tt,1),AllCatchData.pertrial.coord1(tt,2),Catch.ROI1Poly{ii}(:,1),Catch.ROI1Poly{ii}(:,2));
    withinBound2 = inpolygon(AllCatchData.pertrial.coord2(tt,1),AllCatchData.pertrial.coord2(tt,2),Catch.ROI2Poly{ii}(:,1),Catch.ROI2Poly{ii}(:,2));

    if ~isnan(AllCatchData.pertrial.coord1(tt,1)) &&  ~isnan(AllCatchData.pertrial.coord1(tt,2))
        AllCatchData.pertrial.Coord1PolyAcc(tt) =  double(withinBound1);
    else
        AllCatchData.pertrial.Coord1PolyAcc(tt) = NaN;
    end

    if ~isnan(AllCatchData.pertrial.coord2(tt,1)) &&  ~isnan(AllCatchData.pertrial.coord2(tt,2))
        AllCatchData.pertrial.Coord2PolyAcc(tt) =  double(withinBound2);
    else
        AllCatchData.pertrial.Coord2PolyAcc(tt) = NaN;
    end

end


figure
tiledlayout('flow', 'Padding', 'none')
for ag = 1:Ages.NumGroups
    for ii = 1:4

        TwoToneImage = imread(fullfile(main_path, ['/stimuli/CatchTrialStimuli/Twotone/Im_' num2str(40+ii), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);
        TwoToneInsert = insertMatrix(background, TwoToneImage);

        ROI1 = polyshape(Catch.ROI1Poly{ii});
        ROI2 =  polyshape(Catch.ROI2Poly{ii});

        Catch.ROI1Area(ii) = area(ROI1);
        Catch.ROI2Area(ii) = area(ROI2);

        set(gcf,'Color',[1,1,1]);%%%set figure background to white
        x1 = 0:1920; x2 = 0:1080;

        nexttile
        imagesc(TwoToneInsert)
        colormap(gray)
        hold on
        
        %plot ROI 1
        correctID = (AllCatchData.pertrial.ageGroup == ag) & (AllCatchData.pertrial.ImNum == ii) & (AllCatchData.pertrial.Coord1PolyAcc == 1);
        incorrectID = AllCatchData.pertrial.ageGroup == ag & (AllCatchData.pertrial.ImNum == ii) & AllCatchData.pertrial.Coord1PolyAcc == 0;
        scatter(AllCatchData.pertrial.coord1(correctID,1), AllCatchData.pertrial.coord1(correctID,2), sizeMark, [0 0.8 0], 'o', 'filled');
        hold on
        scatter(AllCatchData.pertrial.coord1(incorrectID,1), AllCatchData.pertrial.coord1(incorrectID,2), sizeMark, [0.8 0 0], 'o', 'filled');
        %text(AllCatchData.pertrial.pID(incorrectID),AllCatchData.pertrial.pID(incorrectID), string(Numss))
        hold on

        %plot ROI 2

        correctID = AllCatchData.pertrial.ageGroup == ag & (AllCatchData.pertrial.ImNum == ii) & AllCatchData.pertrial.Coord2PolyAcc == 1;
        incorrectID = AllCatchData.pertrial.ageGroup == ag & (AllCatchData.pertrial.ImNum == ii) & AllCatchData.pertrial.Coord2PolyAcc == 0;
        scatter(AllCatchData.pertrial.coord2(correctID,1), AllCatchData.pertrial.coord2(correctID,2), sizeMark, [0 0.8 0], 'x', 'Linewidth', 2);
        hold on
        scatter(AllCatchData.pertrial.coord2(incorrectID,1), AllCatchData.pertrial.coord2(incorrectID,2), sizeMark, [0.8 0 0], 'x', 'Linewidth', 2);
        %text(AllCatchData.pertrial.pID(incorrectID),AllCatchData.pertrial.pID(incorrectID), string(Numss))
        hold on

        axis off

        title(Ages.GroupNames{ag} + " years")
    end
end

%extract to subs struct

subs.CaCoord1PolyAcc = nan(NumSs,4);
subs.CaCoord2PolyAcc = nan(NumSs,4);

for ss = 1:NumSs
    for ii = 1:4
        ssIdx = ssIDs(ss);
        idx = AllCatchData.pertrial.pID == ssIdx & AllCatchData.pertrial.ImNum == ii;
        if sum(idx) == 1
            subs.CaCoord1PolyAcc(ss,ii) = AllCatchData.pertrial.Coord1PolyAcc(idx);
            subs.CaCoord2PolyAcc(ss,ii) = AllCatchData.pertrial.Coord2PolyAcc(idx);
        elseif sum(idx) > 1 % if there are two of the same images in single particpant
            fprintf('subject %i',ss)
            fprintf('image %i',ii)
            break
        end
    end
end
for ss = 1:NumSs
    ii = 5;
    ssIdx = ssIDs(ss);
    idx = AllCatchData.pertrial.pID == ssIdx & AllCatchData.pertrial.ImNum == ii;
    if sum(idx) == 1
        if isnan(subs.CaCoord1PolyAcc(ss,4))
            subs.CaCoord1PolyAcc(ss,4) = AllCatchData.pertrial.Coord1PolyAcc(idx);
        end

        if isnan(subs.CaCoord2PolyAcc(ss,4))
            subs.CaCoord2PolyAcc(ss,4) = AllCatchData.pertrial.Coord2PolyAcc(idx);
        end

    elseif sum(idx) > 1 % if there are two of the same images in single particpant
        fprintf('subject %i',ss)
        fprintf('image %i',ii)
        break
    end
end

disp('Computing Gray-Scale Polygon Accuracy');

%define variables
subs.GrayScaleCoord1PolyAcc = nan(NumSs,NumIi);
subs.GrayScaleCoord2PolyAcc = nan(NumSs,NumIi);

figure
tiledlayout('flow', 'Padding', 'none', 'TileSpacing', 'tight')
for ag = 1:Ages.NumGroups
    for ii = 1:NumIi
        background = NaN(1080,1920);
        agid = subs.AgeGroup == ag;
        iid = imageIDs(ii);

        % load greyscale
        GrayImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/Greyscale/Im_' num2str(iid), '.jpg']));
        GrayImage = imresize(GrayImage, [680 NaN]);
        GrayImageInsert = insertMatrix(background, GrayImage);

        ROI1 = polyshape(images.ROI1Poly{iid});
        ROI2 = polyshape(images.ROI2Poly{iid});

        %Mean and sigma of age group aiming point for each Grey-scale image
        %ROI 1
        mugs1 = squeeze(mean(subs.GrayScaleCoord1(agid,ii,:),'omitnan'))';
        images.sigmags1x(ii, ag) = cov(subs.GrayScaleCoord1(agid,ii,1));
        images.sigmags1y(ii, ag) = cov(subs.GrayScaleCoord1(agid,ii,2));

        %ROI 2
        mugs2 = squeeze(mean(subs.GrayScaleCoord2(agid,ii,:),'omitnan'))';
        images.sigmags2x(ii, ag) = cov(subs.GrayScaleCoord2(agid,ii,1));
        images.sigmags2y(ii, ag) = cov(subs.GrayScaleCoord2(agid,ii,2));

        % plot Grey-scale distributions
        x1 = 0:1920; x2 = 0:1080;

        nexttile
        imagesc(GrayImageInsert)
        colormap(gray)
        hold on

        plot(ROI1);
        hold on
        plot(ROI2);

        hold on
        axis off
        title(Ages.GroupNames{ag} + " years")

        %Group mean within Polygon?

        withinBound1 = inpolygon(mugs1(1),mugs1(2),images.ROI1Poly{iid}(:,1),images.ROI1Poly{iid}(:,2));
        withinBound2 = inpolygon(mugs2(1),mugs2(2),images.ROI2Poly{iid}(:,1),images.ROI2Poly{iid}(:,2));

        if withinBound1 == 1
            scatter(mugs1(1), mugs1(2), 100, [0 0.8 0], 'o', 'filled');
            hold on
        else
            scatter(mugs1(1), mugs1(2), 100, [0.8 0 0], 'o', 'filled');
            hold on
        end
        if withinBound2
            scatter(mugs2(1), mugs2(2), 100, [0 0.8 0], 'x', 'Linewidth', 2);
            hold on
        else
            scatter(mugs2(1), mugs2(2), 100, [0.8 0 0], 'x','Linewidth', 2);
            hold on
        end

        %Individual within bound?
        for ss = 1 : length(agid)
            if agid(ss)
                withinBound1 = inpolygon(subs.GrayScaleCoord1(ss,ii,1),subs.GrayScaleCoord1(ss,ii,2),images.ROI1Poly{iid}(:,1),images.ROI1Poly{iid}(:,2));
                withinBound2 = inpolygon(subs.GrayScaleCoord2(ss,ii,1),subs.GrayScaleCoord2(ss,ii,2),images.ROI2Poly{iid}(:,1),images.ROI2Poly{iid}(:,2));

                if ~isnan(subs.GrayScaleCoord1(ss,ii,1)) &&  ~isnan(subs.GrayScaleCoord1(ss,ii,2))
                    subs.GrayScaleCoord1PolyAcc(ss,ii) =  double(withinBound1);
                else
                    subs.GrayScaleCoord1PolyAcc(ss,ii) = NaN;
                end

                if ~isnan(subs.GrayScaleCoord2(ss,ii,1)) &&  ~isnan(subs.GrayScaleCoord2(ss,ii,2))
                    subs.GrayScaleCoord2PolyAcc(ss,ii) =  double(withinBound2);
                else
                    subs.GrayScaleCoord2PolyAcc(ss,ii) = NaN;
                end
                %plot ROI 1
                if withinBound1 == 1
                    scatter(subs.GrayScaleCoord1(ss,ii,1), subs.GrayScaleCoord1(ss,ii,2), sizeMark, [0 0.8 0], 'o', 'filled');
                    hold on
                else
                    scatter(subs.GrayScaleCoord1(ss,ii,1), subs.GrayScaleCoord1(ss,ii,2), sizeMark, [0.8 0 0], 'o', 'filled');
                    text(subs.TwoToneCoord1(ss,ii,1),subs.TwoToneCoord1(ss,ii,2), string(ss))
                    hold on
                end
                %plot ROI 2
                if withinBound2
                    scatter(subs.GrayScaleCoord2(ss,ii,1), subs.GrayScaleCoord2(ss,ii,2), sizeMark, [0 0.8 0], 'x', 'Linewidth', 2);
                    hold on
                else
                    scatter(subs.GrayScaleCoord2(ss,ii,1), subs.GrayScaleCoord2(ss,ii,2), sizeMark, [0.8 0 0], 'x', 'Linewidth', 2);
                    text(subs.TwoToneCoord2(ss,ii,1),subs.TwoToneCoord2(ss,ii,2), string(ss))
                    hold on
                end

            end


        end
    end
end

disp('Computing TwoTone Polygon Accuracy');

%define variables
subs.TwoToneCoord1PolyAcc = nan(NumSs,NumIi);
subs.TwoToneCoord2PolyAcc = nan(NumSs,NumIi);


figure
tiledlayout('flow', 'Padding', 'none', 'TileSpacing', 'tight')
for ag = 1:Ages.NumGroups
    for ii = 1:NumIi
        figure
        tiledlayout('flow', 'Padding', 'none', 'TileSpacing', 'tight')

        background = NaN(1080,1920);
        agid = subs.AgeGroup == ag;
        iid = imageIDs(ii);

        TwoToneImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);
        TwoToneInsert = insertMatrix(background, TwoToneImage);

        ROI1 = polyshape(images.ROI1Poly{iid});
        ROI2 = polyshape(images.ROI2Poly{iid});

        %Mean and sigma of age group aiming point for each Two-Tone image
        %ROI 1
        mutt1 = squeeze(mean(subs.TwoToneCoord1(agid,ii,:),'omitnan'))';
        sigmatt1 = cov(subs.TwoToneCoord1(agid,ii,1), subs.TwoToneCoord1(agid,ii,2));
        images.sigmatt1x(ii, ag) = cov(subs.TwoToneCoord1(agid,ii,1));
        images.sigmatt1y(ii, ag) = cov(subs.TwoToneCoord1(agid,ii,2));
        images.ROI1Area(ii) = area(ROI1);
        images.ROI2Area(ii) = area(ROI2);


        %ROI 2
        mutt2 = squeeze(mean(subs.TwoToneCoord2(agid,ii,:),'omitnan'))';
        sigmatt2 = cov(subs.TwoToneCoord2(agid,ii,1), subs.TwoToneCoord2(agid,ii,2));
        images.sigmatt2x(ii, ag) = cov(subs.TwoToneCoord2(agid,ii,1));
        images.sigmatt2y(ii, ag) = cov(subs.TwoToneCoord2(agid,ii,2));


        % plot Two-Tone distributions
        x1 = 0:1920; x2 = 0:1080;

        nexttile
        imagesc(TwoToneInsert)
        colormap(gray)
        hold on

        plot(ROI1);
        hold on
        plot(ROI2);

        hold on
        axis off
        title(Ages.GroupNames{ag} + " years")


        %Group mean within Polygon?

        withinBound1 = inpolygon(mutt1(1),mutt1(2),images.ROI1Poly{iid}(:,1),images.ROI1Poly{iid}(:,2));
        withinBound2 = inpolygon(mutt2(1),mutt2(2),images.ROI2Poly{iid}(:,1),images.ROI2Poly{iid}(:,2));

        if withinBound1 == 1
            scatter(mutt1(1), mutt1(2), 100, [0 0.8 0], 'o', 'filled');
        else
            scatter(mutt1(1), mutt1(2), 100, [0.8 0 0], 'o', 'filled');
        end

        if withinBound2
            scatter(mutt2(1), mutt2(2), 100, [0 0.8 0], 'x', 'Linewidth', 2);
        else
            scatter(mutt2(1), mutt2(2), 100, [0.8 0 0], 'x','Linewidth', 2);
        end

        %Individual within bound?
        for ss = 1 : length(agid)
            if agid(ss)
                withinBound1l = NaN;
                withinBound2l = NaN;
                withinBound1b = NaN;
                withinBound2b = NaN;
                withinBound1 = inpolygon(subs.TwoToneCoord1(ss,ii,1),subs.TwoToneCoord1(ss,ii,2),images.ROI1Poly{iid}(:,1),images.ROI1Poly{iid}(:,2));
                withinBound2 = inpolygon(subs.TwoToneCoord2(ss,ii,1),subs.TwoToneCoord2(ss,ii,2),images.ROI2Poly{iid}(:,1),images.ROI2Poly{iid}(:,2));

                if ~isnan(subs.TwoToneCoord1(ss,ii,1)) &&  ~isnan(subs.TwoToneCoord1(ss,ii,2))
                    subs.TwoToneCoord1PolyAcc(ss,ii) =  double(withinBound1);
                else
                    subs.TwoToneCoord1PolyAcc(ss,ii) = NaN;
                end

                if ~isnan(subs.TwoToneCoord2(ss,ii,1)) &&  ~isnan(subs.TwoToneCoord2(ss,ii,2))
                    subs.TwoToneCoord2PolyAcc(ss,ii) =  double(withinBound2);
                else
                    subs.TwoToneCoord2PolyAcc(ss,ii) = NaN;
                end

                %local polys
                if images.ROI1LocalPoly{iid} ~= images.ROI1Poly{iid}
                    withinBound1l = inpolygon(subs.TwoToneCoord1(ss,ii,1),subs.TwoToneCoord1(ss,ii,2),images.ROI1LocalPoly{iid}(:,1),images.ROI1LocalPoly{iid}(:,2));
                        if subs.TwoToneCoord1PolyAcc(ss,ii) == 0
                            withinBound1b = inpolygon(subs.TwoToneCoord1(ss,ii,1),subs.TwoToneCoord1(ss,ii,2),images.ROI1BigPoly{iid}(:,1),images.ROI1BigPoly{iid}(:,2));
                        else
                            withinBound1b = NaN;
                        end
                    if ~isnan(subs.TwoToneCoord1(ss,ii,1)) &&  ~isnan(subs.TwoToneCoord1(ss,ii,2))
                       
                        subs.TwoToneCoord1LocalPolyAcc(ss,ii) =  double(withinBound1l);
                        subs.TwoToneCoord1BigPolyAcc(ss,ii) = double(withinBound1b);
                        subs.TwoToneCoord1PolyAcc(ss,ii) =  double(withinBound1);
                    else
                        subs.TwoToneCoord1LocalPolyAcc(ss,ii) = NaN;
                        subs.TwoToneCoord1BigPolyAcc(ss,ii) = NaN;
                    end
                else
                    subs.TwoToneCoord1LocalPolyAcc(ss,ii) = NaN;
                    subs.TwoToneCoord1BigPolyAcc(ss,ii) = NaN;
                end

                if images.ROI2LocalPoly{iid} ~= images.ROI2Poly{iid}
                    withinBound2l = inpolygon(subs.TwoToneCoord2(ss,ii,1),subs.TwoToneCoord2(ss,ii,2),images.ROI2LocalPoly{iid}(:,1),images.ROI2LocalPoly{iid}(:,2));
                        if subs.TwoToneCoord2PolyAcc(ss,ii) == 0
                            withinBound2b = inpolygon(subs.TwoToneCoord2(ss,ii,1),subs.TwoToneCoord2(ss,ii,2),images.ROI2BigPoly{iid}(:,1),images.ROI2BigPoly{iid}(:,2));
                        else
                            withinBound2b = NaN;
                        end
                    if ~isnan(subs.TwoToneCoord2(ss,ii,1)) &&  ~isnan(subs.TwoToneCoord2(ss,ii,2))
                        subs.TwoToneCoord2LocalPolyAcc(ss,ii) =  double(withinBound2l);
                        subs.TwoToneCoord2BigPolyAcc(ss,ii) = double(withinBound2b);
                    else
                        subs.TwoToneCoord2LocalPolyAcc(ss,ii) = NaN;
                        subs.TwoToneCoord2BigPolyAcc(ss,ii) = NaN;
                    end
                else
                    subs.TwoToneCoord2LocalPolyAcc(ss,ii) = NaN;
                    subs.TwoToneCoord2BigPolyAcc(ss,ii) = NaN;
                end

                %plot ROI 1
                if withinBound1 == 1
                    scatter(subs.TwoToneCoord1(ss,ii,1), subs.TwoToneCoord1(ss,ii,2), sizeMark, [0 0.8 0], 'o', 'filled');
                    text(subs.TwoToneCoord1(ss,ii,1),subs.TwoToneCoord1(ss,ii,2), string(ss))
                    hold on
                elseif withinBound1l == 1
                    scatter(subs.TwoToneCoord1(ss,ii,1), subs.TwoToneCoord1(ss,ii,2), sizeMark, [0 0 0.8], 'o', 'filled');
                    text(subs.TwoToneCoord1(ss,ii,1),subs.TwoToneCoord1(ss,ii,2), string(ss))
                    hold on
                elseif withinBound1b == 1
                    scatter(subs.TwoToneCoord1(ss,ii,1), subs.TwoToneCoord1(ss,ii,2), sizeMark, [0.8 0 0.8], 'o', 'filled');
                    text(subs.TwoToneCoord1(ss,ii,1),subs.TwoToneCoord1(ss,ii,2), string(ss))
                    hold on
                else
                    scatter(subs.TwoToneCoord1(ss,ii,1), subs.TwoToneCoord1(ss,ii,2), sizeMark, [0.8 0 0], 'o', 'filled');
                    text(subs.TwoToneCoord1(ss,ii,1),subs.TwoToneCoord1(ss,ii,2), string(ss))
                    hold on
                end
                %plot ROI 2
                if withinBound2 == 1
                    scatter(subs.TwoToneCoord2(ss,ii,1), subs.TwoToneCoord2(ss,ii,2), sizeMark, [0 0.8 0], 'x', 'Linewidth', 2);
                    text(subs.TwoToneCoord2(ss,ii,1),subs.TwoToneCoord2(ss,ii,2), string(ss))
                    hold on
                elseif withinBound2l == 1
                    scatter(subs.TwoToneCoord2(ss,ii,1), subs.TwoToneCoord2(ss,ii,2), sizeMark, [0 0 0.8], 'x', 'Linewidth', 2);
                    text(subs.TwoToneCoord2(ss,ii,1),subs.TwoToneCoord2(ss,ii,2), string(ss))
                    hold on
                elseif withinBound2b == 1
                    scatter(subs.TwoToneCoord2(ss,ii,1), subs.TwoToneCoord2(ss,ii,2), sizeMark, [0.8 0 0.8], 'x', 'Linewidth', 2);
                    text(subs.TwoToneCoord2(ss,ii,1),subs.TwoToneCoord2(ss,ii,2), string(ss))
                    hold on
                else
                    scatter(subs.TwoToneCoord2(ss,ii,1), subs.TwoToneCoord2(ss,ii,2), sizeMark, [0.8 0 0], 'x', 'Linewidth', 2);
                    text(subs.TwoToneCoord2(ss,ii,1),subs.TwoToneCoord2(ss,ii,2), string(ss))
                    hold on
                end

            end
        end
        clearvars ROI1l ROI2l
    end
end

%% Local polygon barplots
%Figure 3B

%local errors
clearvars locals incorrects

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

idx1 = ~isnan(mean(subs.TwoToneCoord1LocalPolyAcc, 'omitnan'));
idx2 = ~isnan(mean(subs.TwoToneCoord2LocalPolyAcc, 'omitnan'));

ags = [4 3 2 1];

for ag = 1:Ages.NumGroups
agid = subs.AgeGroup == ags(ag);

    ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:)))| isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0;

    LocalCorrect1 = subs.TwoToneCoord1LocalPolyAcc(agid,idx1);
    LocalCorrect1(ToExclude(:,idx1)) = NaN;
    LocalCorrect1(incorrect(:,idx1)) = NaN;
   
    LocalCorrect2 = subs.TwoToneCoord2LocalPolyAcc(agid,idx2);
    LocalCorrect2(ToExclude(:,idx2)) = NaN;
    LocalCorrect2(incorrect(:,idx2)) = NaN;
   
locals(ag,:) = sum([LocalCorrect1 LocalCorrect2],'omitnan');

    Correct1 = subs.TwoToneCoord1PolyAcc(agid,idx1);
    Correct1(ToExclude(:,idx1)) = NaN;
    Correct1(incorrect(:,idx1)) = NaN;
   
    Correct2 = subs.TwoToneCoord2PolyAcc(agid,idx2);
    Correct2(ToExclude(:,idx2)) = NaN;
    Correct2(incorrect(:,idx2)) = NaN;

incorrects(ag,:) = sum([Correct1==0 Correct2==0],'omitnan');

    BigPoly1 = subs.TwoToneCoord1BigPolyAcc(agid,:);
    BigPoly1(ToExclude) = NaN;
    BigPoly1(incorrect) = NaN;
   
    BigPoly2 = subs.TwoToneCoord2BigPolyAcc(agid,:);
    BigPoly2(ToExclude) = NaN;
    BigPoly2(incorrect) = NaN;


globals(ag,:) = sum([BigPoly1 BigPoly2],'omitnan');

    Poly1 = subs.TwoToneCoord1PolyAcc(agid,:);
    Poly1(ToExclude) = NaN;
    Poly1(incorrect) = NaN;
   
    Poly2 = subs.TwoToneCoord2PolyAcc(agid,:);
    Poly2(ToExclude) = NaN;
    Poly2(incorrect) = NaN;


totalerr(ag,:) = sum([Poly1==0 Poly2==0],'omitnan');
total(ag,:) = sum([~isnan(Poly1) ~isnan(Poly2)],'omitnan');
end

percentlocals = (locals ./ incorrects)*100;
percentglobals = (globals ./ totalerr)*100;
percenttotals = (totalerr ./ total)*100;

%locals plot
b = bar(mean(percentlocals','omitnan'), 'facecolor', 'flat');
hold on
e = errorbar(b.XEndPoints, b.YEndPoints, std(percentlocals','omitnan')/sqrt(length(~isnan(percentlocals)')),...
    'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2, 'LineStyle', 'none');

e.CapSize = 10;
b.CData = [MeanColours{2,1}; MeanColours{2,2}; MeanColours{2,3}; MeanColours{2,4}];

b.LineWidth = 2;
set(gca,'xlim', [0.25 4.75]);
set(gca,'ylim', [0 50]);
set(gca,'xticklabels', {'4-5' '7-9' '10-12' 'adult'}, 'FontSize', 25);
ylabel('Local errors (%)', 'FontSize', 25);
ax = gca;
ax.LineWidth = 2;
box off
ax.XAxis.TickLength = [0 0];


%globals plot
figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

b2 = bar(mean(percentglobals','omitnan'), 'facecolor', 'flat');
hold on
e2 = errorbar(b2.XEndPoints, b2.YEndPoints, std(percentglobals','omitnan')/sqrt(length(~isnan(percentglobals)')),...
    'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2, 'LineStyle', 'none');

e2.CapSize = 10;
b2.CData = [MeanColours{2,1}; MeanColours{2,2}; MeanColours{2,3}; MeanColours{2,4}];

b2.LineWidth = 2;
set(gca,'xlim', [0.25 4.75]);
set(gca,'ylim', [0 50]);
set(gca,'xticklabels', {'4-5' '7-9' '10-12' 'adult'}, 'FontSize', 25);
%set(gca,'ytick', 0:10:100);
ylabel('Global errors (%)', 'FontSize', 25);
ax = gca;
ax.LineWidth = 2;
box off
ax.XAxis.TickLength = [0 0];


%% Image Analysis

if exist('ImageAn.mat', 'file')
    load('ImageAn.mat');
    load('CaImageAn.mat');
    Catch_TTStats = load(fullfile(main_path, '/scripts/Analysis/CatchTwoToneStatsfovG30_fovB30.mat'));
    Catch_GSStats = load(fullfile(main_path, 'scripts/Analysis/CatchGreyScaleStatsfovG30_fovB30.mat'));
    TT_Stats = load(fullfile(main_path, 'scripts/Analysis/TwoToneStatsfovG30_fovB30.mat'));
    GS_Stats = load(fullfile(main_path, 'scripts/Analysis/GreyScaleStatsfovG30_fovB30.mat'));
else
    for ii = 1:NumIi
        iid = imageIDs(ii);
        background = zeros(1080,1920);

        % load two tone
        TwoToneImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);
        TwoToneInsert = insertMatrix(background, TwoToneImage);

        BW1 = edge(TwoToneInsert,'Sobel');
        imshowpair(TwoToneInsert,BW1,'montage')

        %edges of whole image
        imageAn.ID(ii) = iid;
        imageAn.size(ii,:)= size(TwoToneImage);
        imageAn.edges(ii) = (sum(sum(BW1)))./(BGsize(1)*BGsize(2))*100; %calculate edges as percentage of screen size
        imageAn.edges2(ii) = (sum(sum(BW1))./(imageAn.size(ii,1)*imageAn.size(ii,2)))*100; %calculate edges as percentage of image size (edge density)

        % load grey scale
        GreyImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/Greyscale/Im_' num2str(iid), '.jpg']));
        GreyImage = imresize(GreyImage, [680 NaN]);
        GreyImageInsert = insertMatrix(background, GreyImage);

        BW2 = edge(GreyImageInsert,'Sobel');
        imshowpair(GreyImageInsert,BW2,'montage')

        %edges of whole image
        imageAn.GSedges(ii) = (sum(sum(BW2)))./(BGsize(1)*BGsize(2))*100; %calculate edges as percentage of screen size
        imageAn.GSedges2(ii) = (sum(sum(BW2))./(imageAn.size(ii,1)*imageAn.size(ii,2)))*100; %calculate edges as percentage of image size (edge density)
    end

    FourierTT = load(fullfile(main_path,'scripts/Analysis/imstats_tessamooneys/fourierproperties_twotones.mat'));
    WBstatsTT = load(fullfile(main_path, 'scripts/Analysis/imstats_tessamooneys/WBstatistics_defaultsettings_twotones.mat'));
    FourierGS = load(fullfile(main_path, 'scripts/Analysis/imstats_tessamooneys/fourierproperties_grayscale.mat'));
    WBstatsGS = load(fullfile(main_path, 'scripts/Analysis/imstats_tessamooneys/WBstatistics_defaultsettings_grayscale.mat'));

    for ii = 1:NumIi
        for tt = 1:NumIi
            if contains(['Im_' num2str(imageAn.ID(ii))  '.jpg'],FourierTT.names(tt))
                imageAn.TTFourierIntercept(ii) = FourierTT.intercept(tt);
                imageAn.TTFourierSlope(ii) = FourierTT.slope(tt);
            end
            if contains(['Im_' num2str(imageAn.ID(ii))  '.jpg'],FourierGS.names(tt))
                imageAn.GSFourierIntercept(ii) = FourierGS.intercept(tt);
                imageAn.GSFourierSlope(ii) = FourierGS.slope(tt);
            end
            if contains(['Im_' num2str(imageAn.ID(ii))  '.jpg'],WBstatsTT.filenames(tt))
                imageAn.TTLGNBeta(ii) = WBstatsTT.LGNBeta(tt); %Beta parameter varies with the range ofcontrast strengths present in the image (Contrast Energy)
                imageAn.TTLGNGamma(ii) = WBstatsTT.LGNGamma(tt); % Gammaparameter varies with the degree of correlation between contrasts (Spatial Coherance)
                imageAn.TTV1Beta(ii) = WBstatsTT.V1Beta(tt);
                imageAn.TTV1Gamma(ii) = WBstatsTT.V1Gamma(tt);
                imageAn.TTBeta(ii) = WBstatsTT.Beta(tt);
                imageAn.TTGamma(ii) = WBstatsTT.Gamma(tt);
            end
            if contains(['Im_' num2str(imageAn.ID(ii))  '.jpg'],WBstatsGS.filenames(tt))
                imageAn.GSLGNBeta(ii) = WBstatsGS.LGNBeta(tt);
                imageAn.GSLGNGamma(ii) = WBstatsGS.LGNGamma(tt);
                imageAn.GSV1Beta(ii) = WBstatsGS.V1Beta(tt);
                imageAn.GSV1Gamma(ii) = WBstatsGS.V1Gamma(tt);
                imageAn.GSBeta(ii) = WBstatsGS.Beta(tt);
                imageAn.GSGamma(ii) = WBstatsGS.Gamma(tt);
            end
        end
    end
end


%% LOAD Neural network scores from Scores.xls

NetScores = readtable('Scores.xlsx');

% Load 2nd Image set (ImageNet) scores
ImNet.levels = readtable('ImNetScores.xlsx');
ImNet.Ps = readtable('ImNetScores.xlsx', 'Sheet', 'Participants');
ImNet.P1 = readtable('ImNetScores.xlsx', 'Sheet', 'P1');
ImNet.P2 = readtable('ImNetScores.xlsx', 'Sheet', 'P2');
ImNet.P3 = readtable('ImNetScores.xlsx', 'Sheet', 'P3');
ImNet.P4 = readtable('ImNetScores.xlsx', 'Sheet', 'P4');
ImNet.P5 = readtable('ImNetScores.xlsx', 'Sheet', 'P5');
ImNet.P6 = readtable('ImNetScores.xlsx', 'Sheet', 'P6');
ImNet.P7 = readtable('ImNetScores.xlsx', 'Sheet', 'P7');
ImNet.P8 = readtable('ImNetScores.xlsx', 'Sheet', 'P8');
ImNet.P9 = readtable('ImNetScores.xlsx', 'Sheet', 'P9');
ImNet.NasNet = readtable('ImNetScores.xlsx', 'Sheet', 'NasNet');
ImNet.AlexNet= readtable('ImNetScores.xlsx', 'Sheet', 'AlexNet');
ImNet.CorNet = readtable('ImNetScores.xlsx', 'Sheet', 'CorNet-S');


%extract participant (1:9) and Network (10:11) data
ImNet.Ims = ImNet.P1{4:end-1,3}; %order of images presented (all structs in this order)
for i = 1 : length(ImNet.Ims)
    im = ImNet.Ims(i);
    ImNet.smooth(i) = ImNet.levels{im,3};
    ImNet.thresh(i) = ImNet.levels{im,4};
    ImNet.TT1(i,:) = [ImNet.P1{i+3, 6}, ImNet.P2{i+3, 6}, ImNet.P3{i+3, 6}, ImNet.P4{i+3, 6}...
        ImNet.P5{i+3, 6} ImNet.P6{i+3, 6} ImNet.P7{i+3, 6} ImNet.P8{i+3, 6} ImNet.P9{i+3, 6} ...
        ImNet.NasNet{i+3, 6} ImNet.AlexNet{i+3, 6}  ImNet.CorNet{i+3, 6}];
    ImNet.TT2(i,:) = [ImNet.P1{i+3, 8}, ImNet.P2{i+3, 8}, ImNet.P3{i+3, 8}, ImNet.P4{i+3, 8}...
        ImNet.P5{i+3, 8} ImNet.P6{i+3, 8} ImNet.P7{i+3, 8} ImNet.P8{i+3, 8} ImNet.P9{i+3, 8}];
    ImNet.GS(i,:) = [ImNet.P1{i+3, 7}, ImNet.P2{i+3, 7}, ImNet.P3{i+3, 7}, ImNet.P4{i+3, 7}...
        ImNet.P5{i+3, 7} ImNet.P6{i+3, 7} ImNet.P7{i+3, 7} ImNet.P8{i+3, 7} ImNet.P9{i+3, 7} ...
        ImNet.NasNet{i+3, 7} ImNet.AlexNet{i+3, 7} ImNet.CorNet{i+3, 7}];

end

%% Numbers of excluded trials

for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;
    disp(Ages.GroupNames{ag} + " years: " + sum(sum(cell2mat(subs.GreyScaleScored(agid,:)) == 0)) + " excluded trials");
end 

%% Verbal Accuracy Participants (excluding VA-GS inaccuracies) (Figs 1B & S2A)

%change to plot catch or greyscale trials
PlotGS = 1; %Fig1B
PlotCatch = 0; %FigS2A

clearvars CI tmpmean
figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;
    c = linspace(Ages.minmax(ag,1),Ages.minmax(ag,2),sum(agid));

    %to exclude
    ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:)))| isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0;

    %catch trials

        Caincorrect = subs.CaGrayScAcc(agid,:) ~= 1; 

        tmp = subs.CaTwoToneAcc(agid,:);
        tmp(Caincorrect) = NaN;
        Y3 = mean(tmp','omitnan');  
    if PlotCatch
        scatter(subs.ages(agid), Y3,50, [0 0 0], 'LineWidth', 1, 'MarkerEdgeAlpha', 0.3);
        hold on
    end
        %CI95
        tmpMean(ag,2) = mean(Y3,'omitnan'); %#ok<SAGROW> 
        CI(ag,2,:) = bootci(2000,@nanmean,Y3); %#ok<SAGROW>


    %Grey Scale Accuracy

        tmp = cell2mat(subs.GreyScaleScored(agid,:));  
        tmp = double(tmp > 0);
        tmp(ToExclude) = NaN;

        Y1 = mean(tmp','omitnan'); 
    if PlotGS
        %scatter participants
        scatter(subs.ages(agid), Y1,50, 0.5*(subs.ages(agid)), 'filled', 's', 'MarkerFaceAlpha', 0.5); %#ok<UNRCH> 
        hold on
    end
        %CI95
        tmpMean(ag,1) = mean(Y1,'omitnan'); %#ok<SAGROW> 
        CI(ag,1,:) = bootci(2000,@nanmean,Y1); %#ok<SAGROW> 


    %Two Tone Accuracy
    tmp = cell2mat(subs.TwoToneScored(agid,:));
    tmp = double(tmp > 0);
    tmp(ToExclude) = NaN;
    tmp(incorrect) = NaN;


    Y2 = mean(tmp','omitnan');

    %scatter participants
    scatter(subs.ages(agid), Y2, 50, -0.5*subs.ages(agid),'filled', 'MarkerFaceAlpha', .5);
    hold on

    %CI95
    tmpMean(ag,3) = mean(Y2,'omitnan'); %#ok<SAGROW>

    CI(ag,3,:) = bootci(2000,@mean,Y2); %#ok<SAGROW>

    %scatter Catch group means
    if PlotCatch
        errorbar(mean(subs.ages(agid),'omitnan'), tmpMean(ag,2),  tmpMean(ag,2)-CI(ag,2,1), CI(ag,2,2)-tmpMean(ag,2), ...
            'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.5); 
        scatter(mean(subs.ages(agid),'omitnan'), tmpMean(ag,2), 250, [1 1 1],'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);
    end

    %scatter GS group means
    if PlotGS
        errorbar(mean(subs.ages(agid),'omitnan'),tmpMean(ag,1), tmpMean(ag,1)-CI(ag,1,1), CI(ag,1,2)-tmpMean(ag,1), ...
            'Color', [0 0 0],'LineWidth', 1.5); %#ok<UNRCH> 
        scatter(mean(subs.ages(agid),'omitnan'), tmpMean(ag,1),250, 0.5*(mean(subs.ages(agid),'omitnan')), 'filled','s', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);
    end
    hold on

    %scatter TT group means
    errorbar(mean(subs.ages(agid),'omitnan'), tmpMean(ag,3),  tmpMean(ag,3)-CI(ag,3,1), CI(ag,3,2)-tmpMean(ag,3), ...
        'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.5);
    scatter(mean(subs.ages(agid),'omitnan'), tmpMean(ag,3), 250, -0.5*(mean(subs.ages(agid),'omitnan')), 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);
end

set(gca,'LineWidth', 1.5, 'xlim',([0 27.5]),'ylim',([0 1]), 'FontSize', 18);
set(gca,'ytick',([0 0.2 0.4 0.6 0.8 1]), 'yticklabels',{'0', '20', '40', '60', '80', '100'});
set(gca,'xtick',(5:5:25), 'xticklabels', {'5', '10', '15', '20', '25'});
xlabel('Age (years)', 'FontSize', 20);
ylabel('Naive Naming Accuracy (%)', 'FontSize', 20);

if PlotGS
    colormap([spring2;flip(gray(256))]);
else
colormap(spring2(1:end-30,:));
end

%% Neural net TT vs accuracy  (Figs 4A & S2C)

ExclIDs = [4 7 8 11 13 14 16 18 20 21 23 25 28 31 39 40];
idx = mean((ImNet.Ims ~= ExclIDs)') == 1;

%remove koala image  
toDelete = contains(NetScores.Var1,{'Im_8'});
NetScores(toDelete,:) = [];

clearvars  Y Y1G Y1C Y1T
Y1G(1,:) = [mean(NetScores.Alex_GS_Acc), mean(NetScores.CorS_GS_Acc), mean(NetScores.Nas_GS_Acc) 0];

Y1 = NetScores.Alex_TT_Acc(1:4);
Y1(NetScores.Alex_GS_Acc(1:4)==0) = NaN;

Y2 = NetScores.CorS_TT_Acc(1:4);
Y2(NetScores.CorS_GS_Acc(1:4)==0) = NaN;

Y3 = NetScores.Nas_TT_Acc(1:4);
Y3(NetScores.Nas_GS_Acc(1:4)==0) = NaN;

Y1C(1,:) = [mean(Y1, 'omitnan') mean(Y2, 'omitnan') mean(Y3,'omitnan') 0];

Y1 = NetScores.Alex_TT_Acc(5:end);
Y1(NetScores.Alex_GS_Acc(5:end)==0) = NaN;

Y2 = NetScores.CorS_TT_Acc(5:end);
Y2(NetScores.CorS_GS_Acc(5:end)==0) = NaN;

Y3 = NetScores.Nas_TT_Acc(5:end);
Y3(NetScores.Nas_GS_Acc(5:end)==0) = NaN;

Y1T(1,:) = [mean(Y1, 'omitnan') mean(Y2, 'omitnan') mean(Y3,'omitnan') 0];

x = 4;

for ag = [4 3 2 1]
    agid = subs.AgeGroup == ag;
    c = linspace(Ages.minmax(ag,1),Ages.minmax(ag,2),sum(agid));
    x = x + 1;

    %to exclude
    ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:)))| isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0;

    %Grey Scale Accuracy
    tmp = cell2mat(subs.GreyScaleScored(agid,:));
    tmp = double(tmp > 0);
    tmp(ToExclude) = NaN;
    %exclude koala
    tmp(:,6) = NaN;

    Y1 = mean(tmp','omitnan'); 
    Y1G(1, x)= mean(Y1,'omitnan');

    %catch trials
    Caincorrect = subs.CaGrayScAcc(agid,:) ~= 1;

    tmp = subs.CaTwoToneAcc(agid,:);
    tmp(Caincorrect) = NaN;
    Y3 = mean(tmp','omitnan'); 

    Y1C(1, x) = mean(Y3,'omitnan');

    %Two Tone Accuracy
    tmp = cell2mat(subs.TwoToneScored(agid,:));
    tmp = double(tmp > 0);
    tmp(ToExclude) = NaN;
    tmp(incorrect) = NaN;
    %exclude koala
    tmp(:,6) = NaN;

    Y2 = mean(tmp','omitnan'); 
    Y1T(1, x) = mean(Y2,'omitnan');
end

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location
t = tiledlayout(2,2, 'padding', 'normal');

nexttile
b = bar(Y1G, 'FaceColor','flat', 'LineWidth',.8);
hold on

for i = 1:4
    errorbar(b.XEndPoints(9-i), b.YEndPoints(9-i), tmpMean(i,1)-CI(i,1,1), CI(i,1,2)-tmpMean(i,1), ...
        'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.5, 'LineStyle', 'none');
    hold on
end

b.CData = [[1 1 1]; MeanColours{1,1}; MeanColours{1,3}; [0 0 0]; MeanColours{2,1}; MeanColours{2,2}; MeanColours{2,3}; MeanColours{2,4}];

a =[b.XEndPoints(1) b.XEndPoints(2) b.XEndPoints(3) b.XEndPoints(5) b.XEndPoints(6) b.XEndPoints(7) b.XEndPoints(8)];
a = a-0.2;
set(gca,'LineWidth', 1.2, 'xlim',([0 9]),'ylim',([0 1]));
set(gca,'ytick',([0 0.2 0.4 0.6 0.8 1]), 'yticklabels',{'0', '20', '40', '60', '80', '100'});
set(gca,'xtick',a, 'xticklabels', {'AlexNet' 'CorNet' 'NasNet' '4-5' '7-9' '10-12' 'adult' }, 'FontSize', 12);
ylabel('Accuracy (%)', 'FontSize', 14);
h = gca;
h.XAxis.TickLength = [0 0];
box off
xtickangle(45);

nexttile

b = bar(Y1T, 'FaceColor','flat', 'LineWidth',.8);
hold on

for i = 1:4
    errorbar(b.XEndPoints(9-i), b.YEndPoints(9-i), tmpMean(i,3)-CI(i,3,1), CI(i,3,2)-tmpMean(i,3), ...
        'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.5, 'LineStyle', 'none');
    hold on
end

b.CData = [[1 1 1]; MeanColours{1,1}; MeanColours{1,3}; [0 0 0]; MeanColours{2,1}; MeanColours{2,2}; MeanColours{2,3}; MeanColours{2,4}];

a =[b.XEndPoints(1) b.XEndPoints(2) b.XEndPoints(3) b.XEndPoints(5) b.XEndPoints(6) b.XEndPoints(7) b.XEndPoints(8)];
a = a-0.2;
set(gca,'LineWidth', 1.2, 'xlim',([0 9]),'ylim',([0 1]));
set(gca,'ytick',([0 0.2 0.4 0.6 0.8 1]), 'yticklabels',{'0', '20', '40', '60', '80', '100'});
set(gca,'xtick',a, 'xticklabels', {'AlexNet' 'CorNet' 'NasNet' '4-5' '7-9' '10-12' 'adult' }, 'FontSize', 12);
ylabel(["Naive","Accuracy (%)"], 'FontSize', 14);
h = gca;
h.XAxis.TickLength = [0 0];
box off
xtickangle(45);


%DNNs Greyscale
Y(1,1:3) = [mean(ImNet.GS(idx,end-1),'omitnan'), mean(ImNet.GS(idx,end),'omitnan'), mean(ImNet.GS(idx,end-2),'omitnan')];

%Blanks
Y(1,4:7) = 0;
%New Grey Scale Accuracy adults
Y1 = mean(ImNet.GS(idx,1:9)','omitnan'); 

Y(1,8)= mean(Y1,'omitnan');

CI2(1,:)= bootci(2000,@nanmean,Y1); 

%Alexnet
Y1 = ImNet.TT1(idx,end-1);
Y1(ImNet.GS(idx,end-1)==0) = NaN;

%CorNet
Y2 = ImNet.TT1(idx,end);
Y2(ImNet.GS(idx,end)==0) = NaN;

%NasNet
Y3 = ImNet.TT1(idx,end-2);
Y3(ImNet.GS(idx,end-2)==0) = NaN;

Y(2,1:3) = [mean(Y1,'omitnan') mean(Y2, 'omitnan') mean(Y3, 'omitnan')];

%\blanks
Y(2,4:7) = 0;
%Two Tone Accuracy (cued)
Y2 = mean(ImNet.TT2(idx,1:9)','omitnan'); %#ok<UDIM>

Y(2,8) = mean(Y2,'omitnan');

CI2(2,:)= bootci(2000,@nanmean,Y2); %#ok<NANMEAN>

nexttile

b2 = bar(Y(1,:), 'FaceColor','flat', 'EdgeColor',[0 0 0],'LineWidth',1);
hold on


errorbar(b2.XEndPoints(8), b2.YEndPoints(8), b2.YEndPoints(8)-CI2(1,1)', CI2(1,2)'-b2.YEndPoints(8), ...
    'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.5, 'LineStyle', 'none');

hold on


b2.CData = [[1 1 1]; MeanColours{1,1}; MeanColours{1,3}; [1 1 1]; [1 1 1]; [1 1 1]; [1 1 1]; MeanColours{2,4}];

a =[b2.XEndPoints(1) b2.XEndPoints(2) b2.XEndPoints(3) b2.XEndPoints(8)];
a = a-0.2;
set(gca,'LineWidth', 1.2, 'xlim',([0 9]),'ylim',([0 1]));
set(gca,'ytick',([0 0.2 0.4 0.6 0.8 1]), 'yticklabels',{'0', '20', '40', '60', '80', '100'});
set(gca,'xtick',a, 'xticklabels', {'AlexNet' 'CorNet' 'NasNet' 'adult' }, 'FontSize', 12);
ylabel('Accuracy (%)', 'FontSize', 14);
h = gca;
h.XAxis.TickLength = [0 0];
box off
xtickangle(45);

nexttile

b2 = bar(Y(2,:), 'FaceColor','flat', 'EdgeColor',[0 0 0],'LineWidth',1);
hold on

errorbar(b2.XEndPoints(8), b2.YEndPoints(8), b2.YEndPoints(8)-CI2(2,1)', CI2(2,2)'-b2.YEndPoints(8), ...
    'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.5, 'LineStyle', 'none');

hold on

b2.CData = [[1 1 1]; MeanColours{1,1}; MeanColours{1,3}; [1 1 1]; [1 1 1]; [1 1 1]; [1 1 1]; MeanColours{2,4}];

a =[b2.XEndPoints(1) b2.XEndPoints(2) b2.XEndPoints(3) b2.XEndPoints(8)];
a = a-0.2;
set(gca,'LineWidth', 1.2, 'xlim',([0 9]),'ylim',([0 1]));
set(gca,'ytick',([0 0.2 0.4 0.6 0.8 1]), 'yticklabels',{'0', '20', '40', '60', '80', '100'});
set(gca,'xtick',a, 'xticklabels', {'AlexNet' 'CorNet' 'NasNet' 'adult' }, 'FontSize', 12);
ylabel(["Cued/Trained","Accuracy (%)"], 'FontSize', 14);
h = gca;
h.XAxis.TickLength = [0 0];
box off
xtickangle(45);

title(t, 'Greyscale                               Two-tone', 'FontSize',16, 'FontWeight', 'bold')
ylabel(t, 'Trained Image Set           Novel Image Set', 'FontSize',16, 'FontWeight', 'bold')

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location
t = tiledlayout(1,1, 'padding', 'normal');

nexttile

b = bar(Y1C, 'FaceColor','flat', 'LineWidth',1.5);
hold on

for i = 1:4

    errorbar(b.XEndPoints(9-i), b.YEndPoints(9-i), tmpMean(i,1)-CI(i,1,1), CI(i,1,2)-tmpMean(i,1), ...
        'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2, 'LineStyle', 'none');

    hold on

end

b.CData = [[1 1 1]; MeanColours{1,1}; MeanColours{1,3}; [0 0 0]; MeanColours{2,1}; MeanColours{2,2}; MeanColours{2,3}; MeanColours{2,4}];

a =[b.XEndPoints(1) b.XEndPoints(2) b.XEndPoints(3) b.XEndPoints(5) b.XEndPoints(6) b.XEndPoints(7) b.XEndPoints(8)];
a = a-0.2;
set(gca,'LineWidth', 2, 'xlim',([0 9]),'ylim',([0 1]));
set(gca,'ytick',([0 0.2 0.4 0.6 0.8 1]), 'yticklabels',{'0', '20', '40', '60', '80', '100'});
set(gca,'xtick',a, 'xticklabels', {'AlexNet' 'CorNet' 'NasNet' '4-5' '7-9' '10-12' 'adult' }, 'FontSize', 20);
ylabel('Accuracy (%)', 'FontSize', 24);
h = gca;
h.XAxis.TickLength = [0 0];
box off
xtickangle(45);

%% Poly Accuracy data (figs 2B & S2B)
clearvars GSdata TTdata tmpMean CI

PlotGS = 0;
PlotCa = 1;

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    ToExclude = isnan(cell2mat(subs.GreyScaleScored(agid,:))) | isnan(cell2mat(subs.TwoToneScored(agid,:)));

    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,:)) | isnan(subs.TwoToneCoord1PolyAcc(agid,:)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,:)) | isnan(subs.TwoToneCoord2PolyAcc(agid,:));%didn't answer both TT & GS ROI 2

    %catchdata
    Caincorrect = subs.CaGrayScAcc(agid,:) ~= 1;

    if PlotCa

    %ROI 1
    tmp1 = subs.CaCoord1PolyAcc(agid,:);
    tmp1(Caincorrect) = NaN;
    tmp1(mean(ToExclude,2)==1) = NaN;
    tmp1(mean(ToExclude1,2)==1) = NaN;


    %ROI 2
    tmp2 = subs.CaCoord2PolyAcc(agid,:);
    tmp2(Caincorrect) = NaN;
    tmp2(mean(ToExclude,2)==1) = NaN;
    tmp2(mean(ToExclude2,2)==1) = NaN;


    Y3 = mean([tmp1, tmp2]', 'omitnan')'; %#ok<UDIM> %consider ROIs seperately

    scatter(subs.ages(agid), Y3, 50, [0 0 0], 'LineWidth', 1, 'MarkerEdgeAlpha', 0.3)

    hold on

    %CI95

    tmpMean(1) = mean(Y3,'omitnan');

    CI(1,:) = bootci(2000,@nanmean,Y3); %#ok<NANMEAN>
    
    elseif PlotGS

    %Grey Scale Accuracy

    %ROI 1
    tmp1 = subs.GrayScaleCoord1PolyAcc(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN; %answered both TT & GS for each ROI?
    tmp1(incorrect) = NaN; %got GS VA correct?
    %exclude koala
    tmp(:,6) = NaN;

    %ROI 2
    tmp2 = subs.GrayScaleCoord2PolyAcc(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN; %answered both TT & GS for each ROI?
    tmp2(incorrect) = NaN; %got GS VA correct?
    %exclude koala
    tmp(:,6) = NaN;

    Y1 = mean([tmp1, tmp2]', 'omitnan')'; %#ok<UDIM> %consider ROIs seperately

    scatter(subs.ages(agid), Y1,50, 0.5*(subs.ages(agid)), 'filled', 's', 'MarkerFaceAlpha', 0.5); %#ok<UNRCH> 

    hold on

    %CI95
    tmpMean(1) = mean(Y1,'omitnan');

    CI(1,:) = bootci(2000,@nanmean,Y1); %#ok<NANMEAN>
    hold on

    end

    %Two Tone Accuracy
    %ROI 1
    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;
    %exclude koala
    tmp(:,6) = NaN;


    %ROI 2
    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;
    %exclude koala
    tmp(:,6) = NaN;

    Y2 = mean([tmp1, tmp2]', 'omitnan')'; %#ok<UDIM> %consider ROIs seperately

    scatter(subs.ages(agid), Y2, 50, -0.5*subs.ages(agid), 'filled', 'MarkerFaceAlpha', .5);
    hold on

    %CI95

    tmpMean(2) = mean(Y2,'omitnan');

    CI(2,:) = bootci(2000,@nanmean,Y2); %#ok<NANMEAN>


    %catch/GS means
if PlotCa
    h(ag) = errorbar(mean(subs.ages(agid),'omitnan'), tmpMean(1), tmpMean(1)-CI(1,1), CI(1,2)-tmpMean(1), ...
        'Color', [0 0 0],'LineWidth', 1.5);
    scatter(mean(subs.ages(agid),'omitnan'), tmpMean(1),250, [1 1 1], 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);
elseif PlotGS

        h(ag) = errorbar(mean(subs.ages(agid),'omitnan'), tmpMean(1), tmpMean(1)-CI(1,1), CI(1,2)-tmpMean(1), ...
        'Color', [0 0 0],'LineWidth', 1.5);
    scatter(mean(subs.ages(agid),'omitnan'), tmpMean(1),250, 0.5*(mean(subs.ages(agid),'omitnan')), 'filled', 's', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);
end 

    %tt means
    h(ag) = errorbar(mean(subs.ages(agid),'omitnan'), tmpMean(2), tmpMean(2)-CI(2,1), CI(2,2)-tmpMean(2), ...
        'Color', [0 0 0], 'LineWidth', 1.5);
    scatter(mean(subs.ages(agid),'omitnan'), tmpMean(2), 250, -0.5*mean(subs.ages(agid),'omitnan'), 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);
end

set(gca,'LineWidth', 1.5, 'xlim',([0 27.5]),'ylim',([0 1]), 'FontSize', 18);
set(gca,'ytick',([0 0.2 0.4 0.6 0.8 1]), 'yticklabels',{'0', '20', '40', '60', '80', '100'});
set(gca,'xtick', 5:5:25);
xlabel('Age (years)', 'FontSize', 20);
ylabel('Cued Pointing Accuracy (%)', 'FontSize', 20);

if PlotGS
    colormap([spring2;flip(gray(256))]);
else
colormap(spring2(1:end-30,:));
end

%% Poly Accuracy data (Fig 2D)

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report

    TTincorrect = cell2mat(subs.TwoToneScored(agid,:)) == 0; %Twotone report
    TTcorrect = cell2mat(subs.TwoToneScored(agid,:)) >= 1;

    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,:)) | isnan(subs.TwoToneCoord1PolyAcc(agid,:)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,:)) | isnan(subs.TwoToneCoord2PolyAcc(agid,:));%ToExclude = logical(ToExclude + ToExclude2);

    Pincorrect1 = subs.GrayScaleCoord1PolyAcc(agid,:) ~= 1; %GreyScale PolyAcc
    Pincorrect2 = subs.GrayScaleCoord2PolyAcc(agid,:) ~= 1; %GreyScale PolyAcc

    %Two Tone Accuracy (1st TT correct (so remove all TTincorrect), GSVA & GSPoly correct)
    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;
    tmp1(TTincorrect)= NaN;
    tmp1(Pincorrect1) = NaN;

    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;
    tmp2(TTincorrect)= NaN;
    tmp2(Pincorrect2) = NaN;

    %CI95
    Y1 = mean([tmp1, tmp2]', 'omitnan')'; %consider ROIs seperately

    tmpMean1 = mean(Y1,'omitnan');

    scatter(subs.ages(agid), Y1,50, [.5 .5 .5], 'MarkerFaceAlpha', 0.5);
    hold on

    %CI95
    CI1 = bootci(2000,@nanmean,Y1); %#ok<*NANMEAN>

    %Two Tone Accuracy (1st TT incorrect (so remove all TT correct), GSVA & GSPoly correct)
    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;
    tmp1(TTcorrect)= NaN;
    tmp1(Pincorrect1) = NaN;

    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;
    tmp2(TTcorrect)= NaN;
    tmp2(Pincorrect2) = NaN;

    Y2 = mean([tmp1, tmp2]', 'omitnan')'; %#ok<*UDIM> %consider ROIs seperately

    tmpMean = mean(Y2,'omitnan');

    scatter(subs.ages(agid), Y2, 50, 3.14-0.5*subs.ages(agid), 'filled', 'MarkerFaceAlpha', .5);
    hold on

    %CI95
    CI = bootci(2000,@nanmean,Y2);

    h(ag) = errorbar(mean(subs.ages(agid),'omitnan'), tmpMean1, tmpMean1-CI1(1), CI1(2)-tmpMean1, ...
        'Color', [0 0 0], 'LineWidth', 1.5);
    scatter(mean(subs.ages(agid),'omitnan'), tmpMean1, 250, [1 1 1] ,'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);

    h(ag) = errorbar(mean(subs.ages(agid),'omitnan'), tmpMean, tmpMean-CI(1), CI(2)-tmpMean, ...
        'Color', [0 0 0], 'LineWidth', 1.5);
    scatter(mean(subs.ages(agid),'omitnan'), tmpMean, 250, 3.14-0.5*mean(subs.ages(agid),'omitnan'), 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);

    %write out data
    GSdata{ag} = num2str(Y1); %#ok<SAGROW>
    TTdata{ag} = num2str(Y2); %#ok<SAGROW>
end

set(gca,'LineWidth', 1.5, 'xlim',([0 27.5]),'ylim',([0 1]), 'FontSize', 18);
set(gca,'ytick',([0 0.2 0.4 0.6 0.8 1]), 'yticklabels',{'0', '20', '40', '60', '80', '100'});
set(gca,'xtick', 5:5:25);
xlabel('Age (years)', 'FontSize', 20);
ylabel('Cued Pointing Accuracy (%)', 'FontSize', 20);
colormap(spring2(1:end-30,:));

%% Cued Pointing Distance (Fig 2C)

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:))) | isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS

    %Two Tone Accuracy
    tmp1 = subs.TT2GSDist1(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(incorrect) = NaN;

    tmp2 = subs.TT2GSDist1(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(incorrect) = NaN;

    Y = mean([tmp1, tmp2]', 'omitnan')'; %consider ROIs seperately

    tmpMean = mean(Y,'omitnan');
    scatter(subs.ages(agid), Y, 50, 3.14-0.5*subs.ages(agid), 'filled', 'MarkerFaceAlpha', .5);
    hold on

    %CI95

    CI = bootci(2000,@nanmean,Y);

    h(ag) = errorbar(mean(subs.ages(agid),'omitnan'), tmpMean,  tmpMean-CI(1), CI(2)-tmpMean, ...
        'Color', [0 0 0], 'LineWidth', 1.5);
    scatter(mean(subs.ages(agid),'omitnan'), tmpMean, 250, 3.14-0.5*mean(subs.ages(agid),'omitnan'), 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);
end

set(gca,'LineWidth', 1.5, 'xlim',([0 27.5]),'ylim',([0 100]), 'FontSize', 18);
set(gca,'ytick', 0:20:100);
set(gca,'xtick', 5:5:25);
xlabel('Age (years)', 'FontSize', 20);
ylabel('Cued Pointing Distance (mm)', 'FontSize', 20);
colormap(spring2(1:end-30,:));

hold off

%% Cued Pointing Distance fig 2E (perceptual reorganisation)

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report

    TTincorrect = cell2mat(subs.TwoToneScored(agid,:)) == 0; %Twotone report

    ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:))) | isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS

    ToExclude1 = (isnan(subs.GrayScaleCoord1(agid,:,1)) | isnan(subs.GrayScaleCoord1(agid,:,2)))  | (isnan(subs.TwoToneCoord1(agid,:,1)) | isnan(subs.TwoToneCoord1(agid,:,2))); %didn't answer both TT & GS ROI 1
    ToExclude2 = (isnan(subs.GrayScaleCoord2(agid,:,1)) | isnan(subs.GrayScaleCoord2(agid,:,2)))  | (isnan(subs.TwoToneCoord2(agid,:,1)) | isnan(subs.TwoToneCoord2(agid,:,2))); %didn't answer both TT & GS ROI 2

    %Two Tone correct
    tmp1 = subs.TT2GSDist1(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;
    tmp1(TTincorrect) = NaN;

    tmp2 = subs.TT2GSDist2(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;
    tmp2(TTincorrect) = NaN;

    Y1 = mean([tmp1, tmp2]', 'omitnan')'; %consider ROIs seperately

    CI1(:,ag) = bootci(2000,@nanmean,Y1);

    tmpMean1(ag) = mean(Y1,'omitnan');

    scatter(subs.ages(agid), Y1,50,  [.5 .5 .5], 'MarkerFaceAlpha', 0.5);

    hold on

end

clearvars tmp tmp1 tmp2 tmp3

%Twotone Incorrect
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) ~= 1; %GreyScale report

    TTcorrect = cell2mat(subs.TwoToneScored(agid,:)) >= 1;
    ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:))) | isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS


    ToExclude1 = (isnan(subs.GrayScaleCoord1(agid,:,1)) | isnan(subs.GrayScaleCoord1(agid,:,2)))  | (isnan(subs.TwoToneCoord1(agid,:,1)) | isnan(subs.TwoToneCoord1(agid,:,2))); %didn't answer both TT & GS ROI 1
    ToExclude2 = (isnan(subs.GrayScaleCoord2(agid,:,1)) | isnan(subs.GrayScaleCoord2(agid,:,2)))  | (isnan(subs.TwoToneCoord2(agid,:,1)) | isnan(subs.TwoToneCoord2(agid,:,2))); %ToExclude = logical(ToExclude + ToExclude2);


    %Two Tone incorrect
    tmp1 = subs.TT2GSDist1(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;
    tmp1(TTcorrect) = NaN;

    tmp2 = subs.TT2GSDist2(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;
    tmp2(TTcorrect) = NaN;

    Y2 = mean([tmp1, tmp2]', 'omitnan')'; %consider ROIs seperately

    tmpMean2 = mean(Y2,'omitnan');

    CI = bootci(2000,@nanmean,Y2);

    scatter(subs.ages(agid), Y2, 50, 3.14-0.5*subs.ages(agid), 'filled', 'MarkerFaceAlpha', .5);
    hold on



    h(ag) = errorbar(mean(subs.ages(agid),'omitnan'), tmpMean1(ag), tmpMean1(ag)-CI1(1,ag), CI1(2,ag)-tmpMean1(ag), ...
        'Color', [0 0 0], 'LineWidth', 1.5);
    scatter(mean(subs.ages(agid),'omitnan'), tmpMean1(ag),250, [1 1 1], 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);

    hold on

    h(ag) = errorbar(mean(subs.ages(agid),'omitnan'), tmpMean2, tmpMean2-CI(1), CI(2)-tmpMean2, ...
        'Color', [0 0 0], 'LineWidth', 1.5);
    scatter(mean(subs.ages(agid),'omitnan'), tmpMean2,250, 3.14-0.5*mean(subs.ages(agid),'omitnan'), 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);

    %write out data
    TTdata{ag} = num2str(Y2);

end

set(gca,'LineWidth', 1.5, 'xlim',([0 27.5]),'ylim',([0 100]), 'FontSize', 18);
set(gca,'ytick', 0:20:100);
set(gca,'xtick', 5:5:25);
xlabel('Age (years)', 'FontSize', 20);
ylabel('Cued Pointing Distance (mm)', 'FontSize', 20);
colormap(spring2(1:end-30,:));

%% HEATMAP
clearvars heatdata2
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

yvalues2{5-ag} = [Ages.GroupNames{ag} ' Naive']; %#ok<SAGROW> 
heatdata2(5-ag,:) = nan(1,24); %#ok<SAGROW> 

%Catch
Caincorrect = subs.CaGrayScAcc(agid,:) ~= 1;
tmp = subs.CaTwoToneAcc(agid,:);
tmp(Caincorrect) = NaN;
heatdata2(5-ag,1:4) = mean(tmp,'omitnan'); %#ok<SAGROW> 

ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:)))| isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS
incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0;
tmp = cell2mat(subs.TwoToneScored(agid,:));
tmp = double(tmp > 0);
tmp(ToExclude) = NaN;
tmp(incorrect) = NaN;
heatdata2(5-ag,5:end) = mean(tmp, 'omitnan'); %#ok<SAGROW> 
end

%cued
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;
    
    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    Caincorrect = subs.CaGrayScAcc(agid,:) ~= 1;
    
    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,:)) | isnan(subs.TwoToneCoord1PolyAcc(agid,:)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,:)) | isnan(subs.TwoToneCoord2PolyAcc(agid,:));%ToExclude = logical(ToExclude + ToExclude2);
    
    
yvalues2{9-ag} = [Ages.GroupNames{ag} ' Cued'];
    
heatdata2(9-ag,:) = nan(1,24);

%Catch (1st TT incorrect (so remove all TT correct), GSVA & GSPoly correct)
    tmp1 = subs.CaCoord1PolyAcc(agid,:);
    tmp1(Caincorrect) = NaN;
    
    tmp2 = subs.CaCoord1PolyAcc(agid,:);
    tmp2(Caincorrect) = NaN;
    
    Y = (tmp1+tmp2)./2;

heatdata2(9-ag,1:4) = mean(Y,'omitnan');

    %Two Tone Accuracy (1st TT incorrect (so remove all TT correct), GSVA & GSPoly correct)
    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;
    
    
    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;
    
    Y = (tmp1+tmp2)./2;
    
heatdata2(9-ag,5:end) = mean(Y, 'omitnan');
end

%distance^-1
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;
  %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:))) | isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS
    
    %Two Tone Accuracy
    tmp1 = (subs.TT2GSDist1(agid,:)+1).^-1;
    tmp1(ToExclude) = NaN;
    tmp1(incorrect) = NaN;
    
    tmp2 = (subs.TT2GSDist1(agid,:)+1).^-1;
    tmp2(ToExclude) = NaN;
    tmp2(incorrect) = NaN;
    
    Y = (tmp1+tmp2)./2;
    
   yvalues2{13-ag} = [Ages.GroupNames{ag} ' Distance^-1'];
   heatdata2(13-ag,1:4) = NaN;
   heatdata2(13-ag,5:end) = mean(Y, 'omitnan');
    
end

%remove koala
heatdata2 = [heatdata2(:,1:7) heatdata2(:,9:end)];

%Row 5
yvalues2{15} = 'NasNet';
Y1 = NetScores.Nas_TT;
Y1(NetScores.Nas_GS_Acc==0) = NaN;
heatdata2(15,:) = Y1;

%Row 6
yvalues2{14} = 'CorNet-S';
Y1 = NetScores.CorS_TT;
Y1(NetScores.CorS_GS_Acc==0) = NaN;
heatdata2(14,:) = Y1;

%Row 7
yvalues2{13} = 'AlexNet';
Y1 = NetScores.Alex_TT;
Y1(NetScores.Alex_GS_Acc==0) = NaN;
heatdata2(13,:) = Y1;  


yvalues3 = [{'4-5 yrs'} {'7-9 yrs'} {'10-12 yrs'} {'adult'} ...
    {'4-5 yrs'} {'7-9 yrs'} {'10-12 yrs'} {'adult'} {'4-5 yrs'} {'7-9 yrs'} {'10-12 yrs'} {'adult'}...
    {'AlexNet'} {'CorNet'} {'NASNet'}];

figure
set(gcf,'Position',[100,100,650,600]);
map = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];

%Pearsons
[rho, pval] = corr(heatdata2(:,5:23)', heatdata2(:,5:23)', 'rows','pairwise');
sig = pval <= 0.05;
rho2 = rho;
%rho2(~sig)= 0;
h = heatmap(yvalues2, yvalues2, rho2);
h.ColorLimits = [-1 1];

colormap(map)
h.GridVisible = 'off';

%% Spearman rank correlations of image stats to performance
clearvars rho pval ag ranks

%Variables
Names = {'threshold' 'smoothing' 'GS Fourier Intercept' 'GS Fourier Slope' 'TT Fourier Intercept' 'TT Fourier Slope' ...
    'Δ Fourier Intercept' 'Δ Fourier Slope' 'GS LGN Beta' 'GS LGN Gamma' 'TT Beta' 'TT Gamma' 'Δ Beta' 'Δ Gamma'};
X = [imageAn.Cathresh(1:4), imageAn.thresh ;... 'threshold'
    imageAn.Casmooth(1:4) imageAn.smooth ; ... 'smoothing'
    Catch_GSStats.Fi(1:4) GS_Stats.Fi; ...'GS Fourier Intercept'
    Catch_GSStats.Fs(1:4) GS_Stats.Fs;... 'GS Fourier Slope'
    Catch_TTStats.Fi(1:4) TT_Stats.Fi;...'TT Fourier Intercept'
    Catch_TTStats.Fs(1:4) TT_Stats.Fs;...'TT Fourier Slope'
    abs(Catch_TTStats.Fi(1:4)-Catch_GSStats.Fi(1:4)) abs(TT_Stats.Fi-GS_Stats.Fi);... 'Δ Fourier Intercept'
    abs(Catch_TTStats.Fs(1:4)-Catch_GSStats.Fs(1:4)) abs(TT_Stats.Fs-GS_Stats.Fs);... 'Δ Fourier Slope'
    Catch_GSStats.LGNBeta(1:4) GS_Stats.LGNBeta;... 'GS LGN Beta'
    Catch_GSStats.LGNGamma(1:4) GS_Stats.LGNGamma;... 'GS LGN Gamma'
    Catch_TTStats.LGNBeta(1:4) TT_Stats.LGNBeta;...'TT Beta'
    Catch_TTStats.LGNGamma(1:4) TT_Stats.LGNGamma;...'TT Gamma'
    abs(Catch_TTStats.LGNBeta(1:4)-Catch_GSStats.LGNBeta(1:4)) abs(TT_Stats.LGNBeta-GS_Stats.LGNBeta);...'Δ Beta'
    abs(Catch_TTStats.LGNGamma(1:4)-Catch_GSStats.LGNGamma(1:4)) abs(TT_Stats.LGNGamma-GS_Stats.LGNGamma)...'Δ Gamma'
    ];

rho = nan(length(Ages.NumGroups),size(X,1));
pval = nan(length(Ages.NumGroups),size(X,1));

%Naive Accuracy
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    % To exclude = remove trails where either TwoTone or Grayscale is NaN
    ToExclude = [isnan(subs.CaTwoToneAcc(agid,:))| isnan(subs.CaGrayScAcc(agid,:)) isnan(cell2mat(subs.TwoToneScored(agid,:)))| isnan(cell2mat(subs.GreyScaleScored(agid,:)))]; %didn't answer both TT & GS

    % To exclude = remove trails where Grayscale wasn't recognised
    incorrect = [subs.CaGrayScAcc(agid,:) == 0 cell2mat(subs.GreyScaleScored(agid,:)) == 0];

    %Two Tone Accuracy
    tmp = [subs.CaTwoToneAcc(agid,:) cell2mat(subs.TwoToneScored(agid,:))];
    tmp(tmp == 2) = 1;
    tmp(ToExclude) = NaN;
    tmp(incorrect) = NaN;

    Y = mean(tmp,'omitnan');

    %Spearmans rank
    [rho(ag,:), pval(ag,:)] = corr(X',Y','Type', 'Spearman','rows','complete');

end

for ag = 1:4
      
    
    disp("Naive Accuracy: " +  Ages.GroupNames{ag})
  
    Varnames{ag} = matlab.lang.makeValidName(strcat(Ages.GroupNames{ag})); %#ok<SAGROW


    a = [rho(ag,:);pval(ag,:)];
    ranks.(Varnames{ag}) = array2table(a, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(ranks.(Varnames{ag}));
end

clearvars rho pval ag

%Cued
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    caincorrect = subs.CaGrayScAcc(agid,:) ~= 1; %GreyScale report

    ToExclude = isnan(cell2mat(subs.GreyScaleScored(agid,:))) | isnan(cell2mat(subs.TwoToneScored(agid,:)));

    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,:)) | isnan(subs.TwoToneCoord1PolyAcc(agid,:)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,:)) | isnan(subs.TwoToneCoord2PolyAcc(agid,:));%didn't answer both TT & GS ROI 2

    %ROI 1
    cat1 = subs.CaCoord1PolyAcc(agid,:);
    cat1(caincorrect) = NaN;

    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;


    %ROI 2
    cat2 = subs.CaCoord2PolyAcc(agid,:);
    cat2(caincorrect) = NaN;

    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;

    Y = ([cat1 tmp1] + [cat2 tmp2])./2;

    Y = mean(Y,'omitnan');


    %Spearmans rank
    [rho(ag,:), pval(ag,:)] = corr(X',Y','Type', 'Spearman','rows','complete');
end

for ag = 1:4

    disp("Cued Accuracy: " +  Ages.GroupNames{ag})
    Varnames{ag} = matlab.lang.makeValidName(['Cued', strcat(Ages.GroupNames{ag})]);

    a2 = [rho(ag,:);pval(ag,:)];
    ranks.(Varnames{ag}) = array2table(a2, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(ranks.(Varnames{ag}));
end

clearvars rho pval ag

%PR: Correct TT
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    caincorrect = subs.CaGrayScAcc(agid,:) ~= 1; %GreyScale report

    ToExclude = isnan(cell2mat(subs.GreyScaleScored(agid,:))) | isnan(cell2mat(subs.TwoToneScored(agid,:)));

    TTincorrect = cell2mat(subs.TwoToneScored(agid,:)) == 0; %Twotone report
    caTTincorrect = subs.CaTwoToneAcc(agid,:) == 0;
    caTTcorrect = subs.CaTwoToneAcc(agid,:) >=1;

    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,:)) | isnan(subs.TwoToneCoord1PolyAcc(agid,:)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,:)) | isnan(subs.TwoToneCoord2PolyAcc(agid,:));%ToExclude = logical(ToExclude + ToExclude2);

    Pincorrect1 = subs.GrayScaleCoord1PolyAcc(agid,:) ~= 1; %GreyScale PolyAcc
    Pincorrect2 = subs.GrayScaleCoord2PolyAcc(agid,:) ~= 1; %GreyScale PolyAcc

    %ROI 1
    cat1 = subs.CaCoord1PolyAcc(agid,:);
    cat1(caincorrect) = NaN;
    cat1(caTTincorrect) = NaN;

    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;
    tmp1(Pincorrect1) = NaN;
    tmp1(TTincorrect) = NaN;


    %ROI 2
    cat2 = subs.CaCoord2PolyAcc(agid,:);
    cat2(caincorrect) = NaN;
    cat2(caTTincorrect) = NaN;

    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;
    tmp2(Pincorrect2) = NaN;
    tmp2(TTincorrect) = NaN;

    Y = ([cat1 tmp1] + [cat2 tmp2])./2;

    Y = mean(Y,'omitnan');


    %Spearmans rank
    [rho(ag,:), pval(ag,:)] = corr(X',Y','Type', 'Spearman', 'rows','complete');
end
for ag = 1:4

    disp("PR: Correct TT: " +  Ages.GroupNames{ag})

    Varnames{ag} = matlab.lang.makeValidName(['PRCorrect', strcat(Ages.GroupNames{ag})]);

    a2 = [rho(ag,:);pval(ag,:)];
    ranks.(Varnames{ag}) = array2table(a2, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(ranks.(Varnames{ag}));
end

%PR: Inorrect TT
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    caincorrect = subs.CaGrayScAcc(agid,:) ~= 1; %GreyScale report

    ToExclude = isnan(cell2mat(subs.GreyScaleScored(agid,:))) | isnan(cell2mat(subs.TwoToneScored(agid,:)));

    caTTincorrect = subs.CaTwoToneAcc(agid,:) == 0;
    TTcorrect = cell2mat(subs.TwoToneScored(agid,:)) >= 1;
    caTTcorrect = subs.CaTwoToneAcc(agid,:) >=1;

    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,:)) | isnan(subs.TwoToneCoord1PolyAcc(agid,:)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,:)) | isnan(subs.TwoToneCoord2PolyAcc(agid,:));%ToExclude = logical(ToExclude + ToExclude2);

    Pincorrect1 = subs.GrayScaleCoord1PolyAcc(agid,:) ~= 1; %GreyScale PolyAcc
    Pincorrect2 = subs.GrayScaleCoord2PolyAcc(agid,:) ~= 1; %GreyScale PolyAcc

    %ROI 1
    cat1 = subs.CaCoord1PolyAcc(agid,:);
    cat1(caincorrect) = NaN;
    cat1(caTTcorrect) = NaN;

    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;
    tmp1(Pincorrect1) = NaN;
    tmp1(TTcorrect) = NaN;


    %ROI 2
    cat2 = subs.CaCoord2PolyAcc(agid,:);
    cat2(caincorrect) = NaN;
    cat2(caTTcorrect) = NaN;

    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;
    tmp2(Pincorrect2) = NaN;
    tmp2(TTcorrect) = NaN;

    Y = ([cat1 tmp1] + [cat2 tmp2])./2;

    Y = mean(Y,'omitnan');


    %Spearmans rank
    [rho(ag,:), pval(ag,:)] = corr(X',Y','Type', 'Spearman', 'rows','complete');
end
for ag = 1:4

    disp("PR: Incorrect TT: " +  Ages.GroupNames{ag})
    Varnames{ag} = matlab.lang.makeValidName(['PRInorrect', strcat(Ages.GroupNames{ag})]);

    a2 = [rho(ag,:);pval(ag,:)];
    ranks.(Varnames{ag}) = array2table(a2, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(ranks.(Varnames{ag}));
end

clearvars rho pval ag

%Cued Distance (no catch)
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:))) | isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS


    %Two Tone Accuracy
    tmp1 = subs.TT2GSDist1(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(incorrect) = NaN;

    tmp2 = subs.TT2GSDist1(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(incorrect) = NaN;


    Y = (tmp1 + tmp2)./2;

    Y = mean(Y,'omitnan');


    %Spearmans rank
    [rho(ag,:), pval(ag,:)] = corr(X(:,5:end)',Y','Type', 'Spearman','rows','complete');
end
for ag = 1:4
        
    disp("Cued distance: " +  Ages.GroupNames{ag})

    Varnames{ag} = matlab.lang.makeValidName(['Distance', strcat(Ages.GroupNames{ag})]);

    a2 = [rho(ag,:);pval(ag,:)];
    ranks.(Varnames{ag}) = array2table(a2, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(ranks.(Varnames{ag}));
end

clearvars rho pval ag

%Correlate image stats
[rho, pval] = corr(X', X', 'rows','complete');
imagestats.rho = array2table(rho, 'VariableNames', Names, 'RowNames', Names);
imagestats.pval = array2table(pval, 'VariableNames', Names, 'RowNames', Names);


%% Pearsons correlations of image stats to performance
clearvars rho pval ag ranks X Y
%Variables
Names = {'threshold' 'smoothing' 'clutter' 'GS Clutter' 'GS Fourier Intercept' 'GS Fourier Slope' 'TT Fourier Intercept' 'TT Fourier Slope' ...
    'Δ Fourier Intercept' 'Δ Fourier Slope' 'GS LGN Beta' 'GS LGN Gamma' 'TT Beta' 'TT Gamma' 'Δ Beta' 'Δ Gamma' 'Δ Clutter'};
X = [imageAn.Cathresh(1:4) imageAn.thresh ;... 'threshold'
    imageAn.Casmooth(1:4) imageAn.smooth; ... 'smoothing'
    CaimageAn.edges2 imageAn.edges2;... 'clutter'
    CaimageAn.GSedges2 imageAn.GSedges2;... 'GS clutter'
    Catch_GSStats.Fi(1:4) GS_Stats.Fi; ...'GS Fourier Intercept'
    Catch_GSStats.Fs(1:4) GS_Stats.Fs;... 'GS Fourier Slope'
    Catch_TTStats.Fi(1:4) TT_Stats.Fi;...'TT Fourier Intercept'
    Catch_TTStats.Fs(1:4) TT_Stats.Fs;...'TT Fourier Slope'
    abs(Catch_TTStats.Fi(1:4)-Catch_GSStats.Fi(1:4)) abs(TT_Stats.Fi-GS_Stats.Fi);... 'Δ Fourier Intercept'
    abs(Catch_TTStats.Fs(1:4)-Catch_GSStats.Fs(1:4)) abs(TT_Stats.Fs-GS_Stats.Fs);... 'Δ Fourier Slope'
    Catch_GSStats.LGNBeta(1:4) GS_Stats.LGNBeta;... 'GS LGN Beta'
    Catch_GSStats.LGNGamma(1:4) GS_Stats.LGNGamma;... 'GS LGN Gamma'
    Catch_TTStats.LGNBeta(1:4) TT_Stats.LGNBeta;...'TT Beta'
    Catch_TTStats.LGNGamma(1:4) TT_Stats.LGNGamma;...'TT Gamma'
    abs(Catch_TTStats.LGNBeta(1:4)-Catch_GSStats.LGNBeta(1:4)) abs(TT_Stats.LGNBeta-GS_Stats.LGNBeta);...'Δ Beta'
    abs(Catch_TTStats.LGNGamma(1:4)-Catch_GSStats.LGNGamma(1:4)) abs(TT_Stats.LGNGamma-GS_Stats.LGNGamma);...'Δ Gamma'
    CaimageAn.edges2-CaimageAn.GSedges2 imageAn.edges2-imageAn.GSedges2]; %'clutter']

rho = nan(length(Ages.NumGroups),size(X,1));
pval = nan(length(Ages.NumGroups),size(X,1));
%Naive Accuracy
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    % To exclude = remove trails where either TwoTone or Grayscale is NaN
    ToExclude = [isnan(subs.CaTwoToneAcc(agid,:))| isnan(subs.CaGrayScAcc(agid,:)) isnan(cell2mat(subs.TwoToneScored(agid,:)))| isnan(cell2mat(subs.GreyScaleScored(agid,:)))]; %didn't answer both TT & GS
    % To exclude = remove trails where Grayscale wasn't recognised
    incorrect = [subs.CaGrayScAcc(agid,:) == 0 cell2mat(subs.GreyScaleScored(agid,:)) == 0];

    %Two Tone Accuracy
    tmp = [subs.CaTwoToneAcc(agid,:) cell2mat(subs.TwoToneScored(agid,:))];
    tmp(tmp == 2) = 1;
    tmp(ToExclude) = NaN;
    tmp(incorrect) = NaN;

    Y = mean(tmp,'omitnan');

    %Pearsons
    [rho(ag,:), pval(ag,:)] = corr(X',Y','rows','complete');

end

for ag = 1:4
disp("Naive accuracy: " +  Ages.GroupNames{ag})
    Varnames{ag} = matlab.lang.makeValidName(strcat(Ages.GroupNames{ag}));

    a = [rho(ag,:);pval(ag,:)];
    Pearsons.(Varnames{ag}) = array2table(a, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(Pearsons.(Varnames{ag}));
end

clearvars rho pval ag

%Cued
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    caincorrect = subs.CaGrayScAcc(agid,:) ~= 1; %GreyScale report

    ToExclude = isnan(cell2mat(subs.GreyScaleScored(agid,:))) | isnan(cell2mat(subs.TwoToneScored(agid,:)));

    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,:)) | isnan(subs.TwoToneCoord1PolyAcc(agid,:)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,:)) | isnan(subs.TwoToneCoord2PolyAcc(agid,:));%didn't answer both TT & GS ROI 2

    %ROI 1
    cat1 = subs.CaCoord1PolyAcc(agid,:);
    cat1(caincorrect) = NaN;

    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;


    %ROI 2
    cat2 = subs.CaCoord2PolyAcc(agid,:);
    cat2(caincorrect) = NaN;

    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;

    Y = ([cat1 tmp1] + [cat2 tmp2])./2;

    Y = mean(Y,'omitnan');

    %Pearsons
    [rho(ag,:), pval(ag,:)] = corr(X',Y','rows','complete');
end

for ag = 1:4
disp("Cued accuracy: " +  Ages.GroupNames{ag})
    Varnames{ag} = matlab.lang.makeValidName(['Cued', strcat(Ages.GroupNames{ag})]);

    a2 = [rho(ag,:);pval(ag,:)];
    Pearsons.(Varnames{ag}) = array2table(a2, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(Pearsons.(Varnames{ag}));
end

clearvars rho pval ag

%Cued
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    caincorrect = subs.CaGrayScAcc(agid,:) ~= 1; %GreyScale report

    ToExclude = isnan(cell2mat(subs.GreyScaleScored(agid,:))) | isnan(cell2mat(subs.TwoToneScored(agid,:)));

    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,:)) | isnan(subs.TwoToneCoord1PolyAcc(agid,:)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,:)) | isnan(subs.TwoToneCoord2PolyAcc(agid,:));%didn't answer both TT & GS ROI 2

    %ROI 1
    cat1 = subs.CaCoord1PolyAcc(agid,:);
    cat1(caincorrect) = NaN;

    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;

    %ROI 2
    cat2 = subs.CaCoord2PolyAcc(agid,:);
    cat2(caincorrect) = NaN;

    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;

    Y = ([cat1 tmp1] + [cat2 tmp2])./2;

    Y = mean(Y,'omitnan');

    %Pearsons
    [rho(ag,:), pval(ag,:)] = corr(X',Y','rows','complete');
end

for ag = 1:4
disp("Cued accuracy: " +  Ages.GroupNames{ag})
    Varnames{ag} = matlab.lang.makeValidName(['Cued', strcat(Ages.GroupNames{ag})]);

    a2 = [rho(ag,:);pval(ag,:)];
    Pearsons.(Varnames{ag}) = array2table(a2, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(Pearsons.(Varnames{ag}));
end

clearvars rho pval ag

%Nets
clearvars  Y
NetNames = {'Alex', 'CorS', 'Nas'};
%Catch
Y = [NetScores.Alex_TT'; NetScores.CorS_TT'; NetScores.Nas_TT'];
Y([NetScores.Alex_GS_Acc'==0; NetScores.CorS_GS_Acc'==0';NetScores.Nas_GS_Acc'==0]) = NaN;

%remove koala

X = [X(:,1:7) X(:,9:end)];

for ag = 1:3
    %Pearsons
    [rho(ag,:), pval(ag,:)] = corr(X',Y(ag,:)','rows','complete');

    Varnames{ag} = matlab.lang.makeValidName(['Prob', strcat(NetNames{ag})]);

    a2 = [rho(ag,:);pval(ag,:)];
    Pearsons.(Varnames{ag}) = array2table(a2, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(Pearsons.(Varnames{ag}));
end

%PR: Inorrect TT
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    caincorrect = subs.CaGrayScAcc(agid,:) ~= 1; %GreyScale report

    ToExclude = isnan(cell2mat(subs.GreyScaleScored(agid,:))) | isnan(cell2mat(subs.TwoToneScored(agid,:)));

    caTTincorrect = subs.CaTwoToneAcc(agid,:) == 0;
    TTcorrect = cell2mat(subs.TwoToneScored(agid,:)) >= 1;
    caTTcorrect = subs.CaTwoToneAcc(agid,:) >=1;

    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,:)) | isnan(subs.TwoToneCoord1PolyAcc(agid,:)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,:)) | isnan(subs.TwoToneCoord2PolyAcc(agid,:));%ToExclude = logical(ToExclude + ToExclude2);

    Pincorrect1 = subs.GrayScaleCoord1PolyAcc(agid,:) ~= 1; %GreyScale PolyAcc
    Pincorrect2 = subs.GrayScaleCoord2PolyAcc(agid,:) ~= 1; %GreyScale PolyAcc

    %ROI 1
    cat1 = subs.CaCoord1PolyAcc(agid,:);
    cat1(caincorrect) = NaN;
    cat1(caTTcorrect) = NaN;

    tmp1 = subs.TwoToneCoord1PolyAcc(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(ToExclude1) = NaN;
    tmp1(incorrect) = NaN;
    tmp1(Pincorrect1) = NaN;
    tmp1(TTcorrect) = NaN;


    %ROI 2
    cat2 = subs.CaCoord2PolyAcc(agid,:);
    cat2(caincorrect) = NaN;
    cat2(caTTcorrect) = NaN;

    tmp2 = subs.TwoToneCoord2PolyAcc(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(ToExclude2) = NaN;
    tmp2(incorrect) = NaN;
    tmp2(Pincorrect2) = NaN;
    tmp2(TTcorrect) = NaN;

    Y = ([cat1 tmp1] + [cat2 tmp2])./2;

    %remove koala
    Y = [Y(:,1:7) Y(:,9:end)];

    Y = mean(Y,'omitnan');


    %Pearsons
    [rho(ag,:), pval(ag,:)] = corr(X',Y','rows','complete');
end
for ag = 1:4
 disp("PR: Inorrect TT: " +  Ages.GroupNames{ag})
    Varnames{ag} = matlab.lang.makeValidName(['PRInorrect', strcat(Ages.GroupNames{ag})]);

    a2 = [rho(ag,:);pval(ag,:)];
    Pearsons.(Varnames{ag}) = array2table(a2, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(Pearsons.(Varnames{ag}));
end

clearvars rho pval ag

%Cued Distance (no catch)
for ag = 1:Ages.NumGroups
    agid = subs.AgeGroup == ag;

    %to exclude
    incorrect = cell2mat(subs.GreyScaleScored(agid,:)) == 0; %GreyScale report
    ToExclude = isnan(cell2mat(subs.TwoToneScored(agid,:))) | isnan(cell2mat(subs.GreyScaleScored(agid,:))); %didn't answer both TT & GS


    %Two Tone Accuracy
    tmp1 = subs.TT2GSDist1(agid,:);
    tmp1(ToExclude) = NaN;
    tmp1(incorrect) = NaN;

    tmp2 = subs.TT2GSDist1(agid,:);
    tmp2(ToExclude) = NaN;
    tmp2(incorrect) = NaN;


    Y = (tmp1 + tmp2)./2;
    
    %remove koala
    Y = [Y(:,1:7) Y(:,9:end)];

    Y = mean(Y,'omitnan');


    %Pearsons
    [rho(ag,:), pval(ag,:)] = corr(X(:,5:end)',Y','rows','complete');
end
for ag = 1:4
 disp("Cued distance: " +  Ages.GroupNames{ag})
    Varnames{ag} = matlab.lang.makeValidName(['Distance', strcat(Ages.GroupNames{ag})]);

    a2 = [rho(ag,:);pval(ag,:)];
    Pearsons.(Varnames{ag}) = array2table(a2, 'VariableNames', Names, 'RowNames', {'rho' 'p'});

    disp(Pearsons.(Varnames{ag}));
end

clearvars rho pval ag

%Correlate image stats
[rho, pval] = corr(X', X', 'rows','complete');
imagestats.rho = array2table(rho, 'VariableNames', Names, 'RowNames', Names);
imagestats.pval = array2table(pval, 'VariableNames', Names, 'RowNames', Names);

%% Figure 5 Bar plots

%edge density
clearvars CI

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

b = bar([mean([CaimageAn.GSedges2 imageAn.GSedges2]); mean(CaimageAn.edges2); mean(imageAn.edges2)], 'facecolor', 'flat');

CI(1,:) = bootci(2000,@mean,[CaimageAn.GSedges2 imageAn.GSedges2]);
CI(2,:) = bootci(2000,@mean,[CaimageAn.edges2]);
CI(3,:) = bootci(2000,@mean,[imageAn.edges2]);
hold on
e = errorbar(b.XEndPoints, b.YEndPoints, b.YData-CI(:,1)', CI(:,2)'-b.YData, ...
    'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2, 'LineStyle', 'none');
e.CapSize = 10;
b.CData = [0.4470 0.5470 0.7410; MeanColours{2,1}; MeanColours{2,3}];
b.LineWidth = 2;
set(gca,'xlim', [0.25 3.75]);
set(gca,'xticklabels', {'Greyscale' 'Catch' 'Two-tone'}, 'FontSize', 25);
set(gca,'ytick', 0:2:10);
ylabel('Edge Density (%)', 'FontSize', 25);
ax = gca;
ax.LineWidth = 2;
box off

%CE
clearvars CI

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

b = bar([mean([Catch_GSStats.LGNBeta GS_Stats.LGNBeta]); mean(Catch_TTStats.LGNBeta); mean(TT_Stats.LGNBeta)], 'facecolor', 'flat');

CI(1,:) = bootci(2000,@mean,[Catch_GSStats.LGNBeta GS_Stats.LGNBeta]);
CI(2,:) = bootci(2000,@mean,[Catch_TTStats.LGNBeta]);
CI(3,:) = bootci(2000,@mean,[TT_Stats.LGNBeta]);
hold on
e = errorbar(b.XEndPoints, b.YEndPoints, b.YData-CI(:,1)', CI(:,2)'-b.YData, ...
    'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2, 'LineStyle', 'none');
e.CapSize = 10;
b.CData = [0.4470 0.5470 0.7410; MeanColours{2,1}; MeanColours{2,3}];
set(gca,'xlim', [0.25 3.75]);
b.LineWidth = 2;
set(gca,'ylim', [0 0.024])
set(gca,'xticklabels', {'Greyscale' 'Catch' 'Two-tone'}, 'FontSize', 25);
set(gca,'ytick', 0:0.004:0.024);
ylabel('CE', 'FontSize', 25);
ax = gca;
ax.YAxis.Exponent = -3;
ax.LineWidth = 2;
box off

%SC
clearvars CI

figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,500,400]);  %specify figure size and location

b = bar([mean([Catch_GSStats.LGNGamma GS_Stats.LGNGamma]); mean(Catch_TTStats.LGNGamma); mean(TT_Stats.LGNGamma)], 'facecolor', 'flat');

CI(1,:) = bootci(2000,@mean,[Catch_GSStats.LGNGamma GS_Stats.LGNGamma]);
CI(2,:) = bootci(2000,@mean,[Catch_TTStats.LGNGamma]);
CI(3,:) = bootci(2000,@mean,[TT_Stats.LGNGamma]);
hold on
e= errorbar(b.XEndPoints, b.YEndPoints, b.YData-CI(:,1)', CI(:,2)'-b.YData, ...
    'Color', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2, 'LineStyle', 'none');
e.CapSize = 10;
b.CData = [0.4470 0.5470 0.7410; MeanColours{2,1}; MeanColours{2,3}];
b.LineWidth = 2;
set(gca,'ylim', [0 1.2]);
set(gca,'xlim', [0.25 3.75]);
set(gca,'xticklabels', {'Greyscale' 'Catch' 'Two-tone'}, 'FontSize', 25);
set(gca,'ytick', 0:0.4:1.2);
ylabel('SC', 'FontSize', 25);
ax = gca;
ax.LineWidth = 2;
box off

%% SC vs CE (Supplementary figure)

%TwoTone
figure
set(gcf,'Color',[1,1,1])                %set figure background to white
set(gcf,'Position',[100,100,900,900]);  %specify figure size and location

Y = TT_Stats.LGNBeta;

X = TT_Stats.LGNGamma;


for ii = 1:length(X)
    iid = imageIDs(ii);

    if ~isnan(X(ii)) && ~isnan(Y(ii))

        TwoToneImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);

        tmp = size(TwoToneImage);
        ratio = ((tmp(2))*0.95)/(tmp(1)*0.03); %consistnet with ratio of x to y axis

        imagesc([X(ii)-ratio*0.0007, X(ii)+ratio*0.0007], [Y(ii)+0.0007, Y(ii)-0.0007],TwoToneImage)

        hold on

        line([X(ii)-ratio*0.0007,X(ii)+ratio*0.0007,X(ii)+ratio*0.0007,X(ii)-ratio*0.0007,X(ii)-ratio*0.0007],[Y(ii)+0.0007,Y(ii)+0.0007,Y(ii)-0.0007,Y(ii)-0.0007,Y(ii)+0.0007], 'Color', MeanColours{2,3}, 'LineWidth', 2.5)

    end
end

%catch

Y = Catch_TTStats.LGNBeta; %TT Edge Density

X = Catch_TTStats.LGNGamma; %GS Edge Density

for ii = 1:4
    iid = 40+ii;

    if ~isnan(X(ii)) && ~isnan(Y(ii))

        TwoToneImage = imread(fullfile(main_path, ['/stimuli/CatchTrialStimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);

        tmp = size(TwoToneImage);
        ratio = ((tmp(2))*0.95)/(tmp(1)*0.03); %consistnet with ratio of x to y axis

        imagesc([X(ii)-ratio*0.0007, X(ii)+ratio*0.0007], [Y(ii)+0.0007, Y(ii)-0.0007],TwoToneImage)

        hold on

        line([X(ii)-ratio*0.0007,X(ii)+ratio*0.0007,X(ii)+ratio*0.0007,X(ii)-ratio*0.0007,X(ii)-ratio*0.0007],[Y(ii)+0.0007,Y(ii)+0.0007,Y(ii)-0.0007,Y(ii)-0.0007,Y(ii)+0.0007], 'Color', MeanColours{2,1}, 'LineWidth', 2.5)

    end
end

Y = GS_Stats.LGNBeta;

X = GS_Stats.LGNGamma;


for ii = 1:length(X)
    iid = imageIDs(ii);

    if ~isnan(X(ii)) && ~isnan(Y(ii))

        [TwoToneImage, Map, Alpha]  = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/GreyScale/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);

        tmp = size(TwoToneImage);

        ratio = ((tmp(2))*0.95)/(tmp(1)*0.03); %consistnet with ratio of x to y axis

        imagesc([X(ii)-ratio*0.0007, X(ii)+ratio*0.0007], [Y(ii)+0.0007, Y(ii)-0.0007],TwoToneImage)
        colormap(gray)

        hold on
        line([X(ii)-ratio*0.0007,X(ii)+ratio*0.0007,X(ii)+ratio*0.0007,X(ii)-ratio*0.0007,X(ii)-ratio*0.0007],[Y(ii)+0.0007,Y(ii)+0.0007,Y(ii)-0.0007,Y(ii)-0.0007,Y(ii)+0.0007], 'Color', [0.4470 0.5470 0.7410], 'LineWidth', 2.5)

    end
end

%catch

Y = Catch_GSStats.LGNBeta;
X = Catch_GSStats.LGNGamma;

for ii = 1:4
    iid = 40+ii;

    if ~isnan(X(ii)) && ~isnan(Y(ii))

        TwoToneImage = imread(fullfile(main_path, ['/stimuli/CatchTrialStimuli/GreyScale/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);

        tmp = size(TwoToneImage);
        ratio = ((tmp(2))*0.95)/(tmp(1)*0.03); %consistnet with ratio of x to y axis

        imagesc([X(ii)-ratio*0.0007, X(ii)+ratio*0.0007], [Y(ii)+0.0007, Y(ii)-0.0007],TwoToneImage)
        hold on

        line([X(ii)-ratio*0.0007,X(ii)+ratio*0.0007,X(ii)+ratio*0.0007,X(ii)-ratio*0.0007,X(ii)-ratio*0.0007],[Y(ii)+0.0007,Y(ii)+0.0007,Y(ii)-0.0007,Y(ii)-0.0007,Y(ii)+0.0007], 'Color', [0.4470 0.5470 0.7410], 'LineWidth', 2.5)

    end
end

set(gca,'LineWidth', 1.5, 'xlim',([0.35 1.3]),'ylim',([0 0.03]), 'FontSize', 18);
set(gca,'ytick',(0:0.005:0.03))
set(gca,'xtick',(0:0.1:1.3))
xlabel('Spatial Coherence (SC)', 'FontSize', 20);
ylabel('Contrast Energy (CE)', 'FontSize', 20);
set(gca, 'YDir','normal')
box off

%% Perceptual Reorganisation Pointing (fig 3)
%selected images
ims = [1 5 14 40];
groupsinv = [4 3 2 1];
figure
tiledlayout(3,4,'padding', 'none', 'tilespacing', 'none')
background = ones(1080,1920)*255;

%greyscales with all 4-5yo pointing
for ii = 1:length(ims)
    
    iid = ims(ii);
    imsID = imageIDs == iid;

     % load greyscale
        GreyImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/GreyScale/Im_' num2str(iid), '.jpg']));
        GreyImage = imresize(GreyImage, [680 NaN]);
        GreyInsert = insertMatrix(background, GreyImage);
        
        % plot greyscale distributions
        nexttile
            set(gcf,'Color',[1,1,1]);%%%set figure background to white
            x1 = 0:1920; x2 = 0:1080;
            [X1,X2] = meshgrid(x1,x2);
            
            imagesc(GreyInsert)
            colormap('gray')
            hold on
            %4-5 yos 
    agid = subs.AgeGroup == 4;
    incorrect = cell2mat(subs.GreyScaleScored(agid,imsID)) == 0; %GreyScale report
    
    %ROI 1
    tmp1 = [subs.GrayScaleCoord1(agid,imsID,1) subs.GrayScaleCoord1(agid,imsID,2)];
    tmp1(incorrect,:) = NaN;
    
    scatter(tmp1(:,1), tmp1(:,2), 3*sizeMark, MeanColours{2,1}, 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.8, 'LineWidth', 4);
    
    %ROI 2
    tmp1 = [subs.GrayScaleCoord2(agid,imsID,1) subs.GrayScaleCoord2(agid,imsID,2)];
    tmp1(incorrect,:) = NaN;

    scatter(tmp1(:,1), tmp1(:,2), 4*sizeMark, MeanColours{2,1}, 'x', 'MarkerEdgeAlpha', 0.8, 'LineWidth', 3);

    box off
xlim([400 1520])
ylim([150 930])
xticks([]);
yticks([]);
set(gca,'Visible','off')
end

 for ii = 1:length(ims)

        iid = ims(ii);
        imsID = imageIDs == iid;
        
        % load TwoTone
        TwoToneImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);
        TwoToneInsert = insertMatrix(background, TwoToneImage);
        
        % plot Two-tone distributions
            nexttile
            imagesc(TwoToneInsert)
            colormap('gray')
            hold on

    for ag = [4 1]
    agid = subs.AgeGroup == ag;
    incorrect = cell2mat(subs.GreyScaleScored(agid,imsID)) == 0; %GreyScale report
    
    TTincorrect = cell2mat(subs.TwoToneScored(agid,imsID)) == 0; %Twotone report
    TTcorrect = cell2mat(subs.TwoToneScored(agid,imsID)) >= 1;
    
    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,imsID)) | isnan(subs.TwoToneCoord1PolyAcc(agid,imsID)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,imsID)) | isnan(subs.TwoToneCoord2PolyAcc(agid,imsID));%ToExclude = logical(ToExclude + ToExclude2);

    Pincorrect1 = subs.GrayScaleCoord1PolyAcc(agid,imsID) ~= 1; %GreyScale PolyAcc
    Pincorrect2 = subs.GrayScaleCoord2PolyAcc(agid,imsID) ~= 1; %GreyScale PolyAcc
    
    
    %ROI 1
    tmp1 = [subs.TwoToneCoord1(agid,imsID,1) subs.TwoToneCoord1(agid,imsID,2)];
    tmp1(ToExclude1, :) = NaN;
    tmp1(incorrect,:) = NaN;
    tmp1(Pincorrect1,:) = NaN;

    scatter(tmp1(:,1), tmp1(:,2), 3*sizeMark, MeanColours{2, groupsinv(ag)}, 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 3);
    hold on
    end
                    
      
box off
xlim([400 1520])
ylim([150 930])
xticks([])
yticks([])
set(gca,'Visible','off')
 end 

 for ii = 1:length(ims)

        iid = ims(ii);
        imsID = imageIDs == iid;

        % load TwoTone
        TwoToneImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);
        TwoToneInsert = insertMatrix(background, TwoToneImage);
 
        % plot Twotone distributions
            nexttile
            imagesc(TwoToneInsert)
            colormap('gray')
            hold on

    for ag = [4 1]

    agid = subs.AgeGroup == ag;
    incorrect = cell2mat(subs.GreyScaleScored(agid,imsID)) == 0; %GreyScale report
    
    TTincorrect = cell2mat(subs.TwoToneScored(agid,imsID)) == 0; %Twotone report
    TTcorrect = cell2mat(subs.TwoToneScored(agid,imsID)) >= 1;
    
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,imsID)) | isnan(subs.TwoToneCoord2PolyAcc(agid,imsID));%ToExclude = logical(ToExclude + ToExclude2);

    Pincorrect2 = subs.GrayScaleCoord2PolyAcc(agid,imsID) ~= 1; %GreyScale PolyAcc

    %ROI 2
    tmp1 = [subs.TwoToneCoord2(agid,imsID,1) subs.TwoToneCoord2(agid,imsID,2)];
    tmp1(ToExclude2, :) = NaN;
    tmp1(incorrect,:) = NaN;
    tmp1(Pincorrect2,:) = NaN;
    
    scatter(tmp1(:,1), tmp1(:,2), 4*sizeMark, MeanColours{2, groupsinv(ag)}, 'x', 'MarkerEdgeAlpha', 0.8, 'LineWidth', 4);
    hold on
    end             

box off
xlim([400 1520])
ylim([150 930])
xticks([]);
yticks([]);

set(gca,'Visible','off')
 end

%% Perceptual Reorganisation Pointing (sup)

%selected images
ims = imageIDs(mean(imageIDs ~= ims')== 1);
figure
tiledlayout(5,8,'padding', 'tight')
background = zeros(1080,1920)*255;



 for ii = 1:4

        iid = 40+ii;
        
       % load TwoTone
        TwoToneImage = imread(fullfile(main_path, ['/stimuli/CatchTrialStimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);
        TwoToneInsert = insertMatrix(background, TwoToneImage);
        
   %Catch ROI 1     
        % plot Two-tone distributions
            nexttile
            imagesc(TwoToneInsert)
            colormap('gray')
            hold on

    for ag = [4 1]
    agid = subs.AgeGroup == ag;
    incorrect = subs.CaGrayScAcc(agid,ii) == 0;
    
    %ROI 1
    tmp1 = [subs.CaTwoToneCoord1(agid,ii,1) subs.CaTwoToneCoord1(agid,ii,2)];
    tmp1(incorrect,:) = NaN;

    scatter(tmp1(:,1), tmp1(:,2), 3*sizeMark, MeanColours{2, groupsinv(ag)}, 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 3);
    hold on
    end
                    
      
box off
xlim([400 1520])
ylim([150 930])
xticks([])
yticks([])
set(gca,'Visible','off')


 %Catch ROI 2    
        % plot Two-tone distributions
            nexttile
            imagesc(TwoToneInsert)
            colormap('gray')
            hold on

    for ag = [4 1]
    agid = subs.AgeGroup == ag;
    incorrect = subs.CaGrayScAcc(agid,ii) == 0;
    
    %ROI 1
    tmp1 = [subs.CaTwoToneCoord2(agid,ii,1) subs.CaTwoToneCoord2(agid,ii,2)];
    tmp1(incorrect,:) = NaN;

    scatter(tmp1(:,1), tmp1(:,2), 3*sizeMark, MeanColours{2, groupsinv(ag)}, 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 3);
    hold on
    end
                    
      
box off
xlim([400 1520])
ylim([150 930])
xticks([])
yticks([])
set(gca,'Visible','off')
 end 

 for ii = 1:length(ims)
%Twotone ROI 1
        iid = ims(ii);
        imsID = imageIDs == iid;
        
        % load TwoTone
        TwoToneImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);
        TwoToneInsert = insertMatrix(background, TwoToneImage);
        
        % plot Two-tone distributions
            nexttile
            imagesc(TwoToneInsert)
            colormap('gray')
            hold on

    for ag = [4 1]
    agid = subs.AgeGroup == ag;
    incorrect = cell2mat(subs.GreyScaleScored(agid,imsID)) == 0; %GreyScale report
    
    TTincorrect = cell2mat(subs.TwoToneScored(agid,imsID)) == 0; %Twotone report
    TTcorrect = cell2mat(subs.TwoToneScored(agid,imsID)) >= 1;
    
    ToExclude1 = isnan(subs.GrayScaleCoord1PolyAcc(agid,imsID)) | isnan(subs.TwoToneCoord1PolyAcc(agid,imsID)); %didn't answer both TT & GS ROI 1
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,imsID)) | isnan(subs.TwoToneCoord2PolyAcc(agid,imsID));%ToExclude = logical(ToExclude + ToExclude2);

    Pincorrect1 = subs.GrayScaleCoord1PolyAcc(agid,imsID) ~= 1; %GreyScale PolyAcc
    Pincorrect2 = subs.GrayScaleCoord2PolyAcc(agid,imsID) ~= 1; %GreyScale PolyAcc
    
    
    %ROI 1
    tmp1 = [subs.TwoToneCoord1(agid,imsID,1) subs.TwoToneCoord1(agid,imsID,2)];
    tmp1(ToExclude1, :) = NaN;
    tmp1(incorrect,:) = NaN;
    tmp1(Pincorrect1,:) = NaN;

    scatter(tmp1(:,1), tmp1(:,2), 3*sizeMark, MeanColours{2, groupsinv(ag)}, 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 3);
    hold on
    end
                    
      
box off
xlim([400 1520])
ylim([150 930])
xticks([])
yticks([])
set(gca,'Visible','off')

%Twotone ROI 2
        iid = ims(ii);
        imsID = imageIDs == iid;

        % load TwoTone
        TwoToneImage = imread(fullfile(main_path, ['/stimuli/Pilot_Stimuli/TwoTone/Im_' num2str(iid), '.jpg']));
        TwoToneImage = imresize(TwoToneImage, [680 NaN]);
        TwoToneInsert = insertMatrix(background, TwoToneImage);
 
        % plot Twotone distributions
            nexttile
            imagesc(TwoToneInsert)
            colormap('gray')
            hold on

    for ag = [4 1]

    agid = subs.AgeGroup == ag;
    incorrect = cell2mat(subs.GreyScaleScored(agid,imsID)) == 0; %GreyScale report
    
    TTincorrect = cell2mat(subs.TwoToneScored(agid,imsID)) == 0; %Twotone report
    TTcorrect = cell2mat(subs.TwoToneScored(agid,imsID)) >= 1;
    
    ToExclude2 = isnan(subs.GrayScaleCoord2PolyAcc(agid,imsID)) | isnan(subs.TwoToneCoord2PolyAcc(agid,imsID));%ToExclude = logical(ToExclude + ToExclude2);

    Pincorrect2 = subs.GrayScaleCoord2PolyAcc(agid,imsID) ~= 1; %GreyScale PolyAcc

    %ROI 2
    tmp1 = [subs.TwoToneCoord2(agid,imsID,1) subs.TwoToneCoord2(agid,imsID,2)];
    tmp1(ToExclude2, :) = NaN;
    tmp1(incorrect,:) = NaN;
    tmp1(Pincorrect2,:) = NaN;
    
    scatter(tmp1(:,1), tmp1(:,2), 4*sizeMark, MeanColours{2, groupsinv(ag)}, 'x', 'MarkerEdgeAlpha', 0.8, 'LineWidth', 4);
    hold on
    end             

box off
xlim([400 1520])
ylim([150 930])
xticks([]);
yticks([]);

set(gca,'Visible','off')
 end
 


%% Save data
%
nsubs = length(subs.ID);
nIm = NumIi;

variableNames = fieldnames(subs);

T = table(repelem(subs.ID,NumIi), ... % ID
    repelem(subs.ages,NumIi), ... % Age
    repelem(subs.AgeGroup,NumIi), ... % AgeGroup
    reshape(subs.(variableNames{4})', [],1), ... % ImageID
    repmat(imageAn.edges2,1,nsubs)', ... % TT_EdgeDensity
    repmat(imageAn.GSedges2,1,nsubs)', ... % GS_EdgeDensity
    repmat(imageAn.smooth,1,nsubs)', ...% Retreived Smoothing Level
    repmat(imageAn.thresh,1,nsubs)', ... % Retreived Threshold level
    repmat(TT_Stats.Fi,1,nsubs)', ...% TT_Fi (Image analyisis: Two-tone Fourier Intercept)
    repmat(TT_Stats.Fs,1,nsubs)', ... % TT_Fs (Image analyisis: Two-tone Fourier Slope)
    repmat(TT_Stats.LGNBeta,1,nsubs)', ... % TT_LGN_Beta (Image analyisis: Two-tone LGN Beta) 
    repmat(TT_Stats.LGNGamma,1,nsubs)', ... % TT_LGN_Gamma (Image analyisis: Two-tone LGN Gamma) 
    repmat(GS_Stats.Fi,1,nsubs)', ... % GS_Fi (Image analyisis: Grey-Scale Fourier Intercept)
    repmat(GS_Stats.Fs,1,nsubs)', ... % GS_Fs (Image analyisis: Grey-Scale Fourier Slope)
    repmat(GS_Stats.LGNBeta,1,nsubs)', ... % GS_LGN_Beta (Image analyisis: Grey-Scale LGN Beta) 
    repmat(GS_Stats.LGNGamma,1,nsubs)', ... % GS_LGN_Gamma (Image analyisis: Grey-Scale LGN Gamma) 
    zeros(length(repelem(subs.ID,NumIi)),1),... %isCatch
    repmat(images.ROI1Area,1,nsubs)',... % ROI1Area
    repmat(images.ROI2Area,1,nsubs)',... % ROI2Area
    reshape(cell2mat(subs.TwoToneScored)', [],1), ... % Rescored VA_TT
    reshape(cell2mat(subs.GreyScaleScored)', [],1), ... % Rescored VA_GS
    reshape(subs.(variableNames{7})', [],1), ... % can_see_now
    reshape(subs.TT2GSDist1', [],1), ... % Distance_ROI1
    reshape(subs.TT2GSDist2', [],1), ... % Distance_ROI2
    reshape(subs.GrayScaleCoord1PolyAcc', [],1), ... % PA_GS_ROI1
    reshape(subs.GrayScaleCoord2PolyAcc', [],1), ... % PA_GS_ROI2
    reshape(subs.TwoToneCoord1PolyAcc', [],1), ... % PA_TT_ROI1
    reshape(subs.TwoToneCoord2PolyAcc', [],1), ... % PA_TT_ROI2
    'VariableNames', ...
    {'ID', 'Age', 'AgeGroup', 'ImageID', 'TT_EdgeDensity', 'GS_Edge_Density','Smooth', 'Thresh', 'TT_Fi', 'TT_Fs', 'TT_LGN_Beta',...
    'TT_LGN_Gamma','GS_Fi', 'GS_Fs', 'GS_LGN_Beta', 'GS_LGN_Gamma', 'isCatch', 'ROI1Area', 'ROI2Area',...
    'VA_TT', 'VA_GS', 'can_see_now','Distance_ROI1', 'Distance_ROI2','PA_GS_ROI1',...
    'PA_GS_ROI2', 'PA_TT_ROI1', 'PA_TT_ROI2'});

%to save out coordinates, e.g.:
% reshape(subs.(variableNames{8})(:,:,1)', [], 1), reshape(subs.(variableNames{8})(:,:,2)', [], 1), ... % PA_TT_ROI1 [X] + [Y] 

writetable(T,'data_IA4.csv')


nIm = 4;

CatchT = table(repelem(subs.ID,nIm), ... % ID
    repelem(subs.ages,nIm), ... % Age
    repelem(subs.AgeGroup,nIm), ... % AgeGroup
    reshape(subs.CaImNum', [],1), ... % ImageID
    repmat(CaimageAn.edges2,1,nsubs)', ... % TT_EdgeDensity
    repmat(CaimageAn.GSedges2,1,nsubs)', ... % GS_EdgeDensity
    repmat(imageAn.Casmooth(1:4),1,nsubs)', ...% Retreived Smoothing Level
    repmat(imageAn.Cathresh(1:4),1,nsubs)', ... % Retreived Threshold level
    repmat(Catch_TTStats.Fi,1,nsubs)', ...% TT_Fi (Image analyisis: Two-tone Fourier Intercept)
    repmat(Catch_TTStats.Fs,1,nsubs)', ... % TT_Fs (Image analyisis: Two-tone Fourier Slope)
    repmat(Catch_TTStats.LGNBeta,1,nsubs)', ... % TT_LGN_Beta (Image analyisis: Two-tone LGN Beta) 
    repmat(Catch_TTStats.LGNGamma,1,nsubs)', ... % TT_LGN_Gamma (Image analyisis: Two-tone LGN Gamma) 
    repmat(Catch_GSStats.Fi,1,nsubs)', ... % GS_Fi (Image analyisis: Grey-Scale Fourier Intercept)
    repmat(Catch_GSStats.Fs,1,nsubs)', ... % GS_Fs (Image analyisis: Grey-Scale Fourier Slope)
    repmat(Catch_GSStats.LGNBeta,1,nsubs)', ... % GS_LGN_Beta (Image analyisis: Grey-Scale LGN Beta) 
    repmat(Catch_GSStats.LGNGamma,1,nsubs)', ... % GS_LGN_Gamma (Image analyisis: Grey-Scale LGN Gamma) 
    ones(length(repelem(subs.ID,nIm)),1),... %isCatch
    nan(size(repelem(subs.ID,nIm))),... % ROI1Area
    nan(size(repelem(subs.ID,nIm))),... % ROI2Area
    reshape(subs.CaTwoToneAcc', [],1), ... % Rescored VA_TT
    reshape(subs.CaGrayScAcc', [],1), ... % Rescored VA_GS
    reshape(subs.CaCanSeeNow', [],1), ... % can_see_now
    nan(size(repelem(subs.ID,nIm))), ... % Distance_ROI1
    nan(size(repelem(subs.ID,nIm))), ... % Distance_ROI2
    nan(size(repelem(subs.ID,nIm))), ... % PA_GS_ROI1
    nan(size(repelem(subs.ID,nIm))), ... % PA_GS_ROI2
    reshape(subs.CaCoord1PolyAcc', [],1), ... % PA_TT_ROI1
    reshape(subs.CaCoord2PolyAcc', [],1), ... % PA_TT_ROI2
    'VariableNames', ...
    {'ID', 'Age', 'AgeGroup', 'ImageID', 'TT_EdgeDensity', 'GS_Edge_Density','Smooth', 'Thresh', 'TT_Fi', 'TT_Fs', 'TT_LGN_Beta',...
    'TT_LGN_Gamma','GS_Fi', 'GS_Fs', 'GS_LGN_Beta', 'GS_LGN_Gamma', 'isCatch', 'ROI1Area', 'ROI2Area',...
    'VA_TT', 'VA_GS', 'can_see_now','Distance_ROI1', 'Distance_ROI2','PA_GS_ROI1',...
    'PA_GS_ROI2', 'PA_TT_ROI1', 'PA_TT_ROI2'});

writetable(CatchT,'Catch_data.csv')

%write concatenated table
writetable([T ;CatchT], 'data_withCatch.csv')
 
    
