
 clc;
 clear;

% Veri setinin yolu (klasörlerin bulunduğu ana dizin)
datasetPath = 'C:\Users\erenk\Desktop\irisverisetison'; 
outputTrainFolder = 'C:\Users\erenk\Desktop\irisTrainSet';
outputTestFolder = 'C:\Users\erenk\Desktop\irisTestSet';

% Ana çıktı klasörlerini oluştur
if ~exist(outputTrainFolder, 'dir')
    mkdir(outputTrainFolder);
end
if ~exist(outputTestFolder, 'dir')
    mkdir(outputTestFolder);
end

% Ana dizindeki alt klasörleri alın
personFolders = dir(datasetPath);
personFolders = personFolders([personFolders.isdir]); % Sadece klasörleri seç
personFolders = personFolders(~ismember({personFolders.name}, {'.', '..'})); % "." ve ".." hariç

% Her klasörü sırayla işle
for i = 1:length(personFolders)
    personName = personFolders(i).name; % Klasör ismi (kişi adı)
    personPath = fullfile(datasetPath, personName); % Kişinin görüntü klasörü

    % ImageDatastore ile klasörü yükleyin
    personImds = imageDatastore(personPath, ...
        'IncludeSubfolders', false, ...
        'FileExtensions', '.png'); % Sadece .png uzantılı dosyalar

    % Kişinin görüntülerine etiket atayın
    personImds.Labels = repelem(categorical({personName}), numel(personImds.Files), 1);

    % Görüntüleri %70 eğitim, %30 test olarak bölün
    [trainImds, testImds] = splitEachLabel(personImds, 0.7, 'randomized');

    % Eğitim ve test klasörlerine kopyalama
    % Eğitim görüntülerini kopyala
    trainPersonFolder = fullfile(outputTrainFolder, personName);
    if ~exist(trainPersonFolder, 'dir')
        mkdir(trainPersonFolder);
    end
    for j = 1:length(trainImds.Files)
        [~, fileName, ext] = fileparts(trainImds.Files{j});
        copyfile(trainImds.Files{j}, fullfile(trainPersonFolder, [fileName ext]));
    end

    % Test görüntülerini kopyala
    testPersonFolder = fullfile(outputTestFolder, personName);
    if ~exist(testPersonFolder, 'dir')
        mkdir(testPersonFolder);
    end
    for j = 1:length(testImds.Files)
        [~, fileName, ext] = fileparts(testImds.Files{j});
        copyfile(testImds.Files{j}, fullfile(testPersonFolder, [fileName ext]));
    end
end

disp('Veriler başarıyla eğitim ve test kümelerine ayrıldı.');

% Eğitim verilerinin olduğu klasör
trainFolder = 'C:\Users\erenk\Desktop\irisTrainSet'; % Eğitim klasörünüzün yolu
outputProcessedFolder = 'C:\Users\erenk\Desktop\processedTrainImages'; % İşlenmiş görüntüler için yeni klasör

% İşlenmiş görüntüler için çıktı klasörünü oluştur
if ~exist(outputProcessedFolder, 'dir')
    mkdir(outputProcessedFolder);
end

% Eğitim klasöründeki alt klasörleri alın
personFolders = dir(trainFolder);
personFolders = personFolders([personFolders.isdir]); % Sadece klasörleri seç
personFolders = personFolders(~ismember({personFolders.name}, {'.', '..'})); % "." ve ".." hariç

% Her alt klasörü sırayla işle
for i = 1:length(personFolders)
    personName = personFolders(i).name; % Klasör ismi (kişi adı)
    personPath = fullfile(trainFolder, personName); % Kişinin görüntü klasörü

    % Yeni klasörde aynı yapıyı oluştur
    processedPersonFolder = fullfile(outputProcessedFolder, personName);
    if ~exist(processedPersonFolder, 'dir')
        mkdir(processedPersonFolder);
    end

    % Klasördeki görüntüleri alın
    images = dir(fullfile(personPath, '*.png')); % Sadece .png dosyalarını seç

    for j = 1:length(images)
        % Görüntüyü yükle
        imgPath = fullfile(personPath, images(j).name);
        img = imread(imgPath);

        % Gri tonlamaya dönüştür
        if size(img, 3) == 3
            imgGray = rgb2gray(img);
        else
            imgGray = img; % Zaten griyse işlem yapma
        end

        % Gauss filtresi uygula
        imgFiltered = imgaussfilt(imgGray, 2); % Sigma değeri: 2

        % İşlenmiş görüntüyü yeni klasöre kaydet
        processedImgPath = fullfile(processedPersonFolder, images(j).name);
        imwrite(imgFiltered, processedImgPath);
    end
end

disp('Eğitim görüntüleri griye dönüştürüldü, Gauss filtresi uygulandı ve kaydedildi.');

% Canny kenar tespiti
outputCannyFolder = 'C:\Users\erenk\Desktop\cannyEdgeImages'; % Canny sonuçları için yeni klasör

% Canny kenar tespiti sonuçları için çıktı klasörünü oluştur
if ~exist(outputCannyFolder, 'dir')
    mkdir(outputCannyFolder);
end

% İşlenmiş görüntülerdeki alt klasörleri alın
personFolders = dir(outputProcessedFolder);
personFolders = personFolders([personFolders.isdir]);
personFolders = personFolders(~ismember({personFolders.name}, {'.', '..'})); % "." ve ".." hariç

% Her alt klasörü sırayla işle
for i = 1:length(personFolders)
    personName = personFolders(i).name; % Klasör ismi (kişi adı)
    personPath = fullfile(outputProcessedFolder, personName); % İşlenmiş görüntülerin klasörü

    % Yeni klasörde aynı yapıyı oluştur
    cannyPersonFolder = fullfile(outputCannyFolder, personName);
    if ~exist(cannyPersonFolder, 'dir')
        mkdir(cannyPersonFolder);
    end

    % Klasördeki görüntüleri alın
    images = dir(fullfile(personPath, '*.png')); % Sadece .png dosyalarını seç

    for j = 1:length(images)
        % Görüntüyü yükle
        imgPath = fullfile(personPath, images(j).name);
        img = imread(imgPath);

        % Canny kenar tespiti uygula
        edges = edge(img, 'canny');

        % Canny kenar tespit sonuçlarını yeni klasöre kaydet
        cannyImgPath = fullfile(cannyPersonFolder, images(j).name);
        imwrite(edges, cannyImgPath);
    end
end

disp('Canny kenar tespiti başarıyla uygulandı ve sonuçlar kaydedildi.');



% 1. Canny Kenar Tespiti ve Gaussian Filtresi
disp('Canny kenar tespiti ve Gaussian filtresi uygulanıyor...');

% İşlenmiş görüntüler klasöründen gelen çıktı
outputProcessedFolder = 'C:\Users\erenk\Desktop\processedTrainImages'; % İşlenmiş görüntüler
outputCannyFolder = 'C:\Users\erenk\Desktop\cannyEdgeImages'; % Canny sonuçları için klasör
outputFilteredFolder = 'C:\Users\erenk\Desktop\filteredEdgeImages'; % Gaussian sonrası kenar çıktıları

% Gaussian sonrası sonuçları için çıktı klasörünü oluştur
if ~exist(outputFilteredFolder, 'dir')
    mkdir(outputFilteredFolder);
end

% İşlenmiş görüntülerdeki alt klasörleri alın
personFolders = dir(outputProcessedFolder);
personFolders = personFolders([personFolders.isdir]);
personFolders = personFolders(~ismember({personFolders.name}, {'.', '..'}));

% Her alt klasörü sırayla işle
for i = 1:length(personFolders)
    personName = personFolders(i).name; % Klasör ismi (kişi adı)
    personPath = fullfile(outputProcessedFolder, personName); % İşlenmiş görüntülerin klasörü

    % Yeni klasörde aynı yapıyı oluştur
    filteredPersonFolder = fullfile(outputFilteredFolder, personName);
    if ~exist(filteredPersonFolder, 'dir')
        mkdir(filteredPersonFolder);
    end

    % Görüntüleri işleyip kaydet
    images = dir(fullfile(personPath, '*.png'));
    for j = 1:length(images)
        % Görüntüyü yükle
        imgPath = fullfile(personPath, images(j).name);
        img = imread(imgPath);

        % Canny kenar tespiti
        edges = edge(img, 'canny');

        % Gaussian filtresi uygula
        smoothedEdges = imgaussfilt(double(edges), 1); % Sigma değeri: 1
        smoothedEdges = imbinarize(smoothedEdges); % Tekrar ikili hale çevir

        % Gaussian sonrası kenar tespit sonuçlarını kaydet
        filteredImgPath = fullfile(filteredPersonFolder, images(j).name);
        imwrite(smoothedEdges, filteredImgPath);
    end
end

disp('Canny kenar tespiti ve Gaussian filtresi başarıyla tamamlandı.');

% 2. Hough Dönüşümü
disp('Hough dönüşümü uygulanıyor...');

% Hough sonuçları için yeni klasör
outputHoughFolder = 'C:\Users\erenk\Desktop\houghCircleImages';
if ~exist(outputHoughFolder, 'dir')
    mkdir(outputHoughFolder);
end

% Gaussian sonrası kenar klasörlerini alın
personFolders = dir(outputFilteredFolder);
personFolders = personFolders([personFolders.isdir]);
personFolders = personFolders(~ismember({personFolders.name}, {'.', '..'}));

% Her alt klasörü sırayla işle
for i = 1:length(personFolders)
    personName = personFolders(i).name; % Klasör ismi (kişi adı)
    personPath = fullfile(outputFilteredFolder, personName); % Gaussian sonrası kenar klasörü

    % Yeni klasörde aynı yapıyı oluştur
    houghPersonFolder = fullfile(outputHoughFolder, personName);
    if ~exist(houghPersonFolder, 'dir')
        mkdir(houghPersonFolder);
    end

    % Görüntüleri işleyip kaydet
    images = dir(fullfile(personPath, '*.png'));
    for j = 1:length(images)
        % Görüntüyü yükle
        imgPath = fullfile(personPath, images(j).name);
        img = imread(imgPath);

        % Göz bebeği tespiti (küçük yarıçaplar: 15-40 piksel)
        [centersPupil, radiiPupil] = imfindcircles(img, [15 40], ...
            'ObjectPolarity', 'dark', 'Sensitivity', 0.95);

        % Eğer göz bebeği bulunamazsa bir sonraki görüntüye geç
        if isempty(centersPupil)
            disp(['Göz bebeği tespit edilemedi: ', imgPath]);
            continue;
        end

        % Göz bebeğinin merkezi ve çapı
        pupilCenter = centersPupil(1, :); % İlk (en güçlü) daireyi al
        pupilRadius = radiiPupil(1);     % Çap

        % Iris tespiti (büyük yarıçaplar: 50-250 piksel)
        [centersIris, radiiIris] = imfindcircles(img, [50 250], ...
            'ObjectPolarity', 'bright', 'Sensitivity', 0.88);

        % Eğer iris bulunamazsa, göz bebeği merkezine göre manuel iris belirleme
        if isempty(centersIris)
            disp(['İris tespit edilemedi, göz bebeğine göre varsayılan iris çizimi uygulanıyor: ', imgPath]);
            irisCenter = pupilCenter; % Iris merkezi göz bebeğiyle aynı kabul edilir
            irisRadius = 100; % Sabit iris yarıçapı
        else
            % Iris'in merkezi ve çapı
            irisCenter = centersIris(1, :); % İlk (en güçlü) daireyi al
            irisRadius = radiiIris(1);      % Çap
        end

        % Görüntü üzerine göz bebeği ve iris çizimi
        figure('Visible', 'off');
        imshow(img);
        hold on;

        % Göz bebeğini kırmızı daireyle çiz
        viscircles(pupilCenter, pupilRadius, 'Color', 'red', 'LineWidth', 1.5);

        % İrisi yeşil daireyle çiz
        viscircles(irisCenter, irisRadius, 'Color', 'green', 'LineWidth', 1.5);

        hold off;

        % İşlenmiş görüntüyü yeni klasöre kaydet
        houghImgPath = fullfile(houghPersonFolder, images(j).name);
        saveas(gcf, houghImgPath);
        close(gcf);
    end
end

disp('Hough dönüşümü başarıyla tamamlandı ve sonuçlar kaydedildi.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% daugman metodu uygulanacak

% Daugman Metodu uygulanacak görüntüler için klasör

% Hough dönüşüm sonuçlarının bulunduğu klasörün yolu
houghCircleImages = 'C:\Users\erenk\Desktop\houghCircleImages';

% Daugman Metodu ile İris Segmentasyonu
disp('Daugman metodu uygulanıyor...');

% Daugman sonuçları için yeni klasör
outputDaugmanFolder = 'C:\Users\erenk\Desktop\daugmanResults';
if ~exist(outputDaugmanFolder, 'dir')
    mkdir(outputDaugmanFolder);
end

% Hough dönüşüm sonuçları klasörünü kontrol et
personFolders = dir(houghCircleImages);
personFolders = personFolders([personFolders.isdir]);
personFolders = personFolders(~ismember({personFolders.name}, {'.', '..'}));

% Her alt klasörü sırayla işle
for i = 1:length(personFolders)
    personName = personFolders(i).name;
    personPath = fullfile(houghCircleImages, personName);

    % Yeni klasörde aynı yapıyı oluştur
    daugmanPersonFolder = fullfile(outputDaugmanFolder, personName);
    if ~exist(daugmanPersonFolder, 'dir')
        mkdir(daugmanPersonFolder);
    end

    % Görüntüleri işleyip kaydet
    images = dir(fullfile(personPath, '*.png'));
    for j = 1:length(images)
        % Görüntüyü yükle
        imgPath = fullfile(personPath, images(j).name);
        img = imread(imgPath);

        % Daugman metodu ile göz bebeği ve iris segmentasyonu
        [pupilRadius, pupilCenter] = daugmanSegment(img, 15, 40);
        [irisRadius, irisCenter] = daugmanSegment(img, 50, 250);

        % Segmentasyonu görselleştir
        figure('Visible', 'off');
        imshow(img);
        hold on;
        viscircles(pupilCenter, pupilRadius, 'Color', 'red', 'LineWidth', 1.5);
        viscircles(irisCenter, irisRadius, 'Color', 'green', 'LineWidth', 1.5);
        hold off;

        % Segmentasyon sonucunu kaydet
        daugmanImgPath = fullfile(daugmanPersonFolder, images(j).name);
        saveas(gcf, daugmanImgPath);
        close(gcf);
    end
end

disp('Daugman metodu başarıyla uygulandı ve sonuçlar kaydedildi.');

%Daugman segmentasyon fonksiyonu
function [radius, center] = daugmanSegment(img, minRadius, maxRadius)
    [centers, radii] = imfindcircles(img, [minRadius maxRadius], ...
        'ObjectPolarity', 'dark', 'Sensitivity', 0.95);
    if isempty(centers)
        radius = 0;
        center = [0, 0];
    else
        radius = radii(1);
        center = centers(1, :);
    end
end

%%%%%%%

%%gabor dalgacık modeli ile özellik çıkarımı


% Gabor sonuçları için yeni klasör
outputGaborFolder = 'C:\Users\erenk\Desktop\gaborFeatures';
if ~exist(outputGaborFolder, 'dir')
    mkdir(outputGaborFolder);
end

% Daugman metodu sonuç klasöründeki alt klasörleri kontrol et
personFolders = dir(outputDaugmanFolder);
personFolders = personFolders([personFolders.isdir]);
personFolders = personFolders(~ismember({personFolders.name}, {'.', '..'}));

% Her alt klasörü sırayla işle
for i = 1:length(personFolders)
    personName = personFolders(i).name;
    personPath = fullfile(outputDaugmanFolder, personName);

    % Yeni klasörde aynı yapıyı oluştur
    gaborPersonFolder = fullfile(outputGaborFolder, personName);
    if ~exist(gaborPersonFolder, 'dir')
        mkdir(gaborPersonFolder);
    end

    % Görüntüleri işleyip kaydet
    images = dir(fullfile(personPath, '*.png'));
    for j = 1:length(images)
        % Görüntüyü yükle
        imgPath = fullfile(personPath, images(j).name);
        img = imread(imgPath);

        % Gri tonlamaya çevir (Gabor filtresi için)
        if size(img, 3) == 3
            imgGray = rgb2gray(img);
        else
            imgGray = img;
        end

        % Gabor filtrelerini tanımla
        wavelengths = [4, 8, 16]; % Farklı dalga boyları
        orientations = [0, 45, 90, 135]; % Farklı yönelimler

        % Gabor filtrelerini uygula ve sonuçları birleştir
        gaborFeatureMap = zeros(size(imgGray));
        for w = 1:length(wavelengths)
            for o = 1:length(orientations)
                gaborArray = gabor(wavelengths(w), orientations(o));
                gaborMag = imgaborfilt(imgGray, gaborArray);
                gaborFeatureMap = gaborFeatureMap + gaborMag; % Özellik haritalarını birleştir
            end
        end

        % Özellik haritasını normalize et
        gaborFeatureMap = mat2gray(gaborFeatureMap);

        % Sonuçları kaydet
        gaborImgPath = fullfile(gaborPersonFolder, images(j).name);
        imwrite(gaborFeatureMap, gaborImgPath);
    end
end

disp('Gabor dalgacıkları başarıyla uygulandı ve sonuçlar kaydedildi.');


%%%% vektör analiz makine öğrenmesi

%Eğitim ve test veri klasörlerinin yolu
trainFolder = 'C:\Users\erenk\Desktop\irisTrainSet';
testFolder = 'C:\Users\erenk\Desktop\irisTestSet';
%Eğitim verisini dengeli hale getirme

%Eğitim verilerinden özellik çıkarma
disp('Eğitim verilerinden özellik çıkarılıyor...');
[trainFeatures, trainLabels] = extractFeaturesFromFolder(trainFolder);
disp('Eğitim verilerinden özellik çıkarma tamamlandı.');

%Test verilerinden özellik çıkarma
disp('Test verilerinden özellik çıkarılıyor...');
[testFeatures, testLabels] = extractFeaturesFromFolder(testFolder);
disp('Test verilerinden özellik çıkarma tamamlandı.');

%Özellik matrisini ve etiketleri doğru formata dönüştür
trainFeatures = double(trainFeatures); % Özellikler double olmalı
testFeatures = double(testFeatures);

%Çoklu sınıflar için SVM Modelini Eğit
disp('SVM modeli (fitcecoc) eğitiliyor...');
svmModel = fitcecoc(trainFeatures, trainLabels, 'Learners', 'linear', 'Coding', 'onevsall');
disp('SVM modeli başarıyla eğitildi.');

%SVM Modelini Test Et
disp('SVM modeli test ediliyor...');
predictions = predict(svmModel, testFeatures);

%Doğruluk oranını hesapla
accuracy = sum(predictions == testLabels) / numel(testLabels);
fprintf('SVM model doğruluğu: %.2f%%\n', accuracy * 100);

%Hata Matrisi
figure;
confusionchart(testLabels, predictions);
title('Hata Matrisi (Confusion Matrix)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% buraya kadar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hata yok
%Yeni Görüntüler İçin Tahmin (Örnek Kullanım)
newImagePath = 'C:\Users\erenk\Desktop\newIris.png'; % Yeni iris görüntüsünün yolu
newImage = imread(newImagePath);
newFeatureVector = double(newImage(:)'); % Özellik vektörüne dönüştür
predictedLabel = predict(svmModel, newFeatureVector);
disp(['Tahmin edilen kişi: ', char(predictedLabel)]);

%Özellik çıkarma fonksiyonu
function [features, labels] = extractFeaturesFromFolder(folderPath)
    features = [];
    labels = [];

    % Klasördeki alt klasörleri alın
    personFolders = dir(folderPath);
    personFolders = personFolders([personFolders.isdir]);
    personFolders = personFolders(~ismember({personFolders.name}, {'.', '..'}));

    % Her alt klasörü sırayla işle
    for i = 1:length(personFolders)
        personName = personFolders(i).name;
        personPath = fullfile(folderPath, personName);
        images = dir(fullfile(personPath, '*.png'));

        for j = 1:length(images)
            % Görüntüyü yükle
            imgPath = fullfile(personPath, images(j).name);
            img = imread(imgPath);

            % Görüntüyü vektöre dönüştür
            featureVector = img(:)'; % Satır vektörüne dönüştür
            features = [features; featureVector];

            % Etiketi ekle
            labels = [labels; {personName}];
        end
    end

    % Etiketleri kategorik hale getir
    labels = categorical(labels);
end

%%%%%% duyarlılık kavramlarının hesaplanmaası

% Performans metriklerini hesaplamak için fonksiyon
function evaluateModel(testLabels, predictions)
    % Hata matrisini oluştur
    confusionMat = confusionmat(testLabels, predictions); % Hata matrisi
    numClasses = size(confusionMat, 1); % Sınıf sayısı

    % Toplam örnek sayısı
    totalSamples = sum(confusionMat(:));

    % Başlangıç değerleri
    TP = zeros(numClasses, 1); % Doğru Pozitif
    FN = zeros(numClasses, 1); % Yanlış Negatif
    FP = zeros(numClasses, 1); % Yanlış Pozitif
    TN = zeros(numClasses, 1); % Doğru Negatif

    % Her sınıf için metrikleri hesapla
    for i = 1:numClasses
        TP(i) = confusionMat(i, i); % Doğru Pozitif
        FN(i) = sum(confusionMat(i, :)) - TP(i); % Yanlış Negatif
        FP(i) = sum(confusionMat(:, i)) - TP(i); % Yanlış Pozitif
        TN(i) = totalSamples - (TP(i) + FP(i) + FN(i)); % Doğru Negatif
    end

    % Doğruluk (Accuracy)
    accuracy = sum(TP) / totalSamples;
    fprintf('Doğruluk (Accuracy): %.2f%%\n', accuracy * 100);

    % Her sınıf için duyarlılık (Recall/Sensitivity)
    recall = TP ./ (TP + FN);
    recall(isnan(recall)) = 0; % NaN durumlarını 0 yap
    fprintf('Duyarlılık (Recall/Sensitivity) (Ortalama): %.2f%%\n', mean(recall) * 100);

    % Her sınıf için kesinlik (Precision)
    precision = TP ./ (TP + FP);
    precision(isnan(precision)) = 0; % NaN durumlarını 0 yap
    fprintf('Kesinlik (Precision) (Ortalama): %.2f%%\n', mean(precision) * 100);

    % Her sınıf için özgüllük (Specificity)
    specificity = TN ./ (TN + FP);
    specificity(isnan(specificity)) = 0; % NaN durumlarını 0 yap
    fprintf('Özgüllük (Specificity) (Ortalama): %.2f%%\n', mean(specificity) * 100);

    % Sınıf bazında sonuçlar
    fprintf('\nSınıf Bazında Metrikler:\n');
    for i = 1:numClasses
        fprintf('Sınıf %d - Recall: %.2f%%, Precision: %.2f%%, Specificity: %.2f%%\n', ...
            i, recall(i) * 100, precision(i) * 100, specificity(i) * 100);
    end
end

% testLabels ve predictions değişkenlerini kontrol edin
disp('Test Etiketleri (Gerçek Değerler):');
disp(testLabels);

disp('Model Tahminleri:');
disp(predictions);

% Performans metriklerini hesapla
disp('Performans metrikleri hesaplanıyor...');
evaluateModel(testLabels, predictions);