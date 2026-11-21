%% ================================================================
% Far-field polarization + BIC topological charge
%
% Outputs:
%   1. Polarization ellipse map
%   2. S3/S0 map
%   3. Ellipticity-angle chi map
%   4. BIC candidate from minimum S0
%   5. Polarization topological charge q
%
% Color convention:
%   S3 > 0  -> Blue
%   S3 < 0  -> Red
%   chi > 0 -> Blue
%   chi < 0 -> Red
%% ================================================================

clear;
clc;
close all;

%% ================================================================
% 1. Read data
%% ================================================================

filename = fullfile(fileparts(mfilename('fullpath')), 'beita1 x.txt');

data = readmatrix(filename, 'FileType','text', 'CommentStyle','%');
assert(size(data,2)>=5,'Expected kx, ky, frequency (THz), cx, cy.');

c = 3e8;
a = 320e-9;

unit = 2*pi*c/a;

kkx  = data(:,1);
kky  = data(:,2);

freq = data(:,3)*1e12;

cx = data(:,4);
cy = data(:,5);


%% ================================================================
% 2. Normalize frequency and select target band
%% ================================================================

freqNorm = 2*pi*real(freq)/unit;

freqMin = 0.5;
freqMax = 0.66;

mask_f = ...
    (freqNorm >= freqMin) & ...
    (freqNorm <= freqMax);

kkx = kkx(mask_f);
kky = kky(mask_f);

cx = cx(mask_f);
cy = cy(mask_f);

freqNorm = freqNorm(mask_f);


%% ================================================================
% 3. Far-field Stokes parameters
%% ================================================================

S0 = abs(cx).^2 + abs(cy).^2;

S1 = abs(cx).^2 - abs(cy).^2;

S2 = 2*real(cx .* conj(cy));

S3 = -2*imag(cx .* conj(cy));


%% ================================================================
% 4. Remove invalid points
%% ================================================================

valid = isfinite(kkx) & isfinite(kky) & isfinite(S0);
assert(any(valid),'No finite samples in the frequency window.');
% Retain radiation zeros for candidate detection and winding interpolation.
allK = [kkx(valid), kky(valid)];
allS0 = S0(valid);
allS1 = S1(valid);
allS2 = S2(valid);
assert(size(unique(allK,'rows'),1)==size(allK,1), ...
    'Repeated k points: select one mode per k before interpolation.');
valid = valid & S0 > 0;
points = [kkx(valid),kky(valid)];
assert(size(points,1)>=3 && rank(points-mean(points,1))==2, ...
    'Polarization maps require noncollinear radiating samples.');

kkx = kkx(valid);
kky = kky(valid);

cx = cx(valid);
cy = cy(valid);

S0 = S0(valid);
S1 = S1(valid);
S2 = S2(valid);
S3 = S3(valid);

freqNorm = freqNorm(valid);


%% ================================================================
% 5. Normalized Stokes parameters
%% ================================================================

S1_norm = S1 ./ S0;
S2_norm = S2 ./ S0;
S3_norm = S3 ./ S0;


%% ================================================================
% 6. Polarization parameters
%% ================================================================

% Polarization azimuth
phi = 0.5 * atan2(S2,S1);

% Signed ellipticity angle
chi = 0.5 * atan2(S3, sqrt(S1.^2 + S2.^2));

% Ellipticity ratio
ellipticity = tan(abs(chi));


%% ================================================================
% 7. Print basic polarization information
%% ================================================================

fprintf('\n================================================\n');
fprintf('        FAR-FIELD POLARIZATION ANALYSIS\n');
fprintf('================================================\n');

fprintf('Number of selected points = %d\n',length(kkx));

[maxS3,idxMax] = max(abs(S3_norm));

fprintf('\nMaximum |S3/S0| point:\n');
fprintf('kx = %.8f\n',kkx(idxMax));
fprintf('ky = %.8f\n',kky(idxMax));
fprintf('S3/S0 = %.8f\n',S3_norm(idxMax));
fprintf('chi = %.6f deg\n',chi(idxMax)*180/pi);
fprintf('ellipticity b/a = %.6f\n',ellipticity(idxMax));


%% ================================================================
% 8. Define common color map
%
% Need:
%   Negative -> Red
%   Zero     -> Middle
%   Positive -> Blue
%
% flipud(jet) gives this orientation.
%% ================================================================

commonCmap = flipud(jet(256));


%% ================================================================
% 9. Polarization ellipse map
%% ================================================================

figure;
hold on;

theta = linspace(0,2*pi,100);

scaleFac = 0.0008;

for cnt = 1:length(kkx)

    phiTemp = phi(cnt);

    chiTemp = abs(chi(cnt));

    kxTemp = kkx(cnt);
    kyTemp = kky(cnt);

    % ------------------------------------------------------------
    % Major axis
    % ------------------------------------------------------------

    majorAxis = scaleFac;

    % Minor / major = tan(|chi|)
    minorAxis = scaleFac * tan(chiTemp);

    % Avoid numerical problems
    minorAxis = max(minorAxis,0);

    % Ellipse before rotation
    x = majorAxis*cos(theta);
    y = minorAxis*sin(theta);

    % Rotate ellipse
    x1 = x*cos(phiTemp) - y*sin(phiTemp) + kxTemp;
    y1 = x*sin(phiTemp) + y*cos(phiTemp) + kyTemp;


    % ------------------------------------------------------------
    % Color classification
    % ------------------------------------------------------------

    if abs(S3_norm(cnt)) < 0.01

        % Approximately linear polarization
        dx = scaleFac*cos(phiTemp);
        dy = scaleFac*sin(phiTemp);

        plot( ...
            [kxTemp-dx,kxTemp+dx], ...
            [kyTemp-dy,kyTemp+dy], ...
            'k', ...
            'LineWidth',2);


    elseif S3_norm(cnt) > 0

        % Positive S3 -> BLUE
        plot( ...
            x1,y1, ...
            'b', ...
            'LineWidth',2);


    else

        % Negative S3 -> RED
        plot( ...
            x1,y1, ...
            'r', ...
            'LineWidth',2);

    end

end

hold off;

axis equal;
grid on;

xlabel('k_x (2*pi/a)');
ylabel('k_y (2*pi/a)');

title('Far-field polarization ellipses');
% ===== 坐标范围 =====

% ===== 让坐标区域成为正方形 =====
axis equal;

% ===== 坐标轴位置：正方形 =====


%% ================================================================
% 10. Generate interpolation grid
%% ================================================================

xmin = min(kkx);
xmax = max(kkx);

ymin = min(kky);
ymax = max(kky);

Nx = 500;
Ny = 500;

[Xg,Yg] = meshgrid( ...
    linspace(xmin,xmax,Nx), ...
    linspace(ymin,ymax,Ny));


%% ================================================================
% 11. Interpolate S3/S0
%% ================================================================

F_S3 = scatteredInterpolant( ...
    kkx,kky,S3_norm, ...
    'natural','none');

S3_grid = F_S3(Xg,Yg);


%% ================================================================
% 12. S3/S0 map
%% ================================================================

figure;

imagesc( ...
    Xg(1,:), ...
    Yg(:,1), ...
    S3_grid);

axis xy equal tight;

% IMPORTANT:
% Same color direction as polarization ellipse map
% Negative -> Red
% Positive -> Blue
colormap(commonCmap);

caxis([-1 1]);

cb = colorbar;

cb.Label.String = 'S_3/S_0';

xlabel('k_x (2\pi/a)');
ylabel('k_y (2\pi/a)');

title('Far-field circular polarization S_3/S_0');


%% ================================================================
% 13. Mark near-circular polarization regions
%% ================================================================

thresholdCircular = 0.99;

maskCircular = ...
    abs(S3_grid) >= thresholdCircular;

[row,col] = find(maskCircular);

hold on;

scatter( ...
    Xg(1,col), ...
    Yg(row,1), ...
    8, ...
    'w', ...
    'filled');

hold off;


%% ================================================================
% 14. Interpolate ellipticity angle chi
%% ================================================================

F_chi = scatteredInterpolant( ...
    kkx,kky,chi, ...
    'natural','none');

chi_grid = F_chi(Xg,Yg);


%% ================================================================
% 15. Signed ellipticity-angle map
%% ================================================================

figure;

imagesc( ...
    Xg(1,:), ...
    Yg(:,1), ...
    chi_grid);

axis xy equal tight;

% Same color direction
% chi < 0 -> Red
% chi > 0 -> Blue
colormap(commonCmap);

caxis([-pi/4 pi/4]);

cb = colorbar;

cb.Ticks = ...
    [-pi/4 -pi/8 0 pi/8 pi/4];

cb.TickLabels = ...
    {'-\pi/4','-\pi/8','0','\pi/8','\pi/4'};

cb.Label.String = '\chi';

xlabel('k_x (2\pi/a)');
ylabel('k_y (2\pi/a)');

title('Far-field ellipticity angle \chi');


%% ================================================================
% 16. Find BIC candidate from minimum far-field radiation
%% ================================================================

% BIC condition:
%
% S0 = |cx|^2 + |cy|^2 -> 0

[minS0,idxBIC] = min(allS0);

kxBIC = allK(idxBIC,1);
kyBIC = allK(idxBIC,2);

fprintf('\n================================================\n');
fprintf('               BIC CANDIDATE\n');
fprintf('================================================\n');

fprintf('kx_BIC = %.10f\n',kxBIC);
fprintf('ky_BIC = %.10f\n',kyBIC);

fprintf('S0_min = %.8e\n',minS0);

%% ================================================================
% 17. Plot BIC candidate
%% ================================================================

figure;

scatter(kkx,kky,15,'k','filled');

hold on;

plot( ...
    kxBIC,kyBIC, ...
    'ro', ...
    'MarkerSize',10, ...
    'LineWidth',2);

xlabel('k_x (2\pi/a)');
ylabel('k_y (2\pi/a)');

axis equal;
grid on;

title('BIC candidate: minimum S_0');

legend( ...
    'Far-field data', ...
    'BIC candidate');


%% ================================================================
% 18. Interpolants for S1 and S2
%% ================================================================

F_S1 = scatteredInterpolant( ...
    allK(:,1),allK(:,2),allS1, ...
    'natural','none');

F_S2 = scatteredInterpolant( ...
    allK(:,1),allK(:,2),allS2, ...
    'natural','none');


%% ================================================================
% 19. Topological charge around BIC
%% ================================================================

radii = [0.002 0.003 0.004 0.005 0.006];

Nloop = 600;

thetaLoop = linspace(0,2*pi,Nloop+1);

qList = nan(size(radii));

fprintf('\n================================================\n');
fprintf('          BIC TOPOLOGICAL CHARGE\n');
fprintf('================================================\n');

for ir = 1:length(radii)

    radius = radii(ir);

    % ------------------------------------------------------------
    % Integration circle
    % ------------------------------------------------------------

    kxLoop = ...
        kxBIC + radius*cos(thetaLoop);

    kyLoop = ...
        kyBIC + radius*sin(thetaLoop);


    % ------------------------------------------------------------
    % Interpolation
    % ------------------------------------------------------------

    S1_loop = F_S1(kxLoop,kyLoop);
    S2_loop = F_S2(kxLoop,kyLoop);


    % ------------------------------------------------------------
    % Valid points
    % ------------------------------------------------------------

    validLoop = ...
        ~isnan(S1_loop) & ...
        ~isnan(S2_loop);

    nValid = sum(validLoop);

    fprintf('\nr = %.5f\n',radius);

    fprintf('valid points = %d / %d\n', ...
        nValid,length(thetaLoop));


    % Require a closed loop with defined polarization azimuth.
    if ~all(validLoop) || any(hypot(S1_loop,S2_loop) <= 64*eps(max(hypot(S1_loop,S2_loop))))

        fprintf( ...
            'This radius is NOT sufficiently covered by data.\n');

        continue;

    end


    % ------------------------------------------------------------
    % Use only valid points
    % ------------------------------------------------------------

    thetaValid = thetaLoop(validLoop);

    S1Valid = S1_loop(validLoop);
    S2Valid = S2_loop(validLoop);


    % ------------------------------------------------------------
    % Polarization azimuth
    % ------------------------------------------------------------

    phiLoop = ...
        0.5*atan2(S2Valid,S1Valid);


    % ------------------------------------------------------------
    % Unwrap 2phi
    %
    % Polarization angle is defined modulo pi
    % ------------------------------------------------------------

    phase2 = unwrap(2*phiLoop);

    phiUnwrapped = phase2/2;


    % ------------------------------------------------------------
    % Winding
    % ------------------------------------------------------------

    deltaPhi = ...
        phiUnwrapped(end) - phiUnwrapped(1);

    q = deltaPhi/(2*pi);

    qList(ir) = q;

    fprintf( ...
        'Delta phi = %.8f rad\n',deltaPhi);

    fprintf( ...
        'q = %.8f\n',q);

end


%% ================================================================
% 20. Stable topological charge
%% ================================================================

validQ = qList(~isnan(qList));

if isempty(validQ)

    warning( ...
        ['No valid integration radius. ' ...
         'The BIC is probably too close to the data boundary, ' ...
         'or the k-space sampling does not cover the integration loops.']);

    qBIC = NaN;

else

    % Median of valid results
    if numel(validQ)>=2 && max(validQ)-min(validQ)<1e-5
        qBIC = median(validQ);
    else
        warning('Insufficient agreement across integration radii. Inspect qList.');
        qBIC = NaN;
    end

end


fprintf('\n-----------------------------------------------\n');

fprintf( ...
    'Estimated BIC topological charge q = %.6f\n', ...
    qBIC);

fprintf('-----------------------------------------------\n');


%% ================================================================
% 21. Plot integration loop
%% ================================================================

if ~isnan(qBIC)

    figure;

    scatter( ...
        kkx,kky, ...
        12, ...
        'k', ...
        'filled');

    hold on;


    % ------------------------------------------------------------
    % Use first valid radius
    % ------------------------------------------------------------

    validIndex = find(~isnan(qList),1,'first');

    radiusPlot = radii(validIndex);


    thetaPlot = linspace(0,2*pi,601);


    kxLoopPlot = ...
        kxBIC + radiusPlot*cos(thetaPlot);

    kyLoopPlot = ...
        kyBIC + radiusPlot*sin(thetaPlot);


    % ------------------------------------------------------------
    % Integration loop
    % ------------------------------------------------------------

    plot( ...
        kxLoopPlot, ...
        kyLoopPlot, ...
        'r-', ...
        'LineWidth',2);


    % ------------------------------------------------------------
    % BIC point
    % ------------------------------------------------------------

    plot( ...
        kxBIC,kyBIC, ...
        'bo', ...
        'MarkerFaceColor','b', ...
        'MarkerSize',8);


    axis equal;
    grid on;

    xlabel('k_x (2\pi/a)');
    ylabel('k_y (2\pi/a)');

    title(sprintf( ...
        'BIC polarization loop, q = %.3f', ...
        qBIC));

    legend( ...
        'Far-field data', ...
        'Integration loop', ...
        'BIC candidate');

    hold off;

end


%% ================================================================
% 22. Plot polarization-angle winding
%% ================================================================

if isfinite(qBIC)

    % ------------------------------------------------------------
    % Interpolate S1 and S2 on selected loop
    % ------------------------------------------------------------

    S1_loop = ...
        F_S1(kxLoopPlot,kyLoopPlot);

    S2_loop = ...
        F_S2(kxLoopPlot,kyLoopPlot);


    % ------------------------------------------------------------
    % Polarization angle
    % ------------------------------------------------------------

    phiLoop = ...
        0.5*atan2(S2_loop,S1_loop);


    % ------------------------------------------------------------
    % IMPORTANT:
    % Polarization angle is modulo pi.
    % Therefore unwrap 2*phi instead of phi directly.
    % ------------------------------------------------------------

    phase2 = unwrap(2*phiLoop);

    phiUnwrapped = phase2/2;


    %% ------------------------------------------------------------
    % Plot
    %% ------------------------------------------------------------

    figure;

    plot( ...
        thetaPlot, ...
        phiUnwrapped, ...
        'LineWidth',2);

    grid on;
    box on;


    %% ------------------------------------------------------------
    % X-axis
    %% ------------------------------------------------------------

    xlim([0 2*pi]);

    xticks(0:pi/2:2*pi);

    xticklabels({ ...
        '0', ...
        '\pi/2', ...
        '\pi', ...
        '3\pi/2', ...
        '2\pi'});

    xlabel('\theta');


    %% ------------------------------------------------------------
    % Y-axis
    %% ------------------------------------------------------------

    phiMin = min(phiUnwrapped);
    phiMax = max(phiUnwrapped);

    yMin = ...
        floor(phiMin/(pi/2))*(pi/2);

    yMax = ...
        ceil(phiMax/(pi/2))*(pi/2);

    if yMin == yMax, yMax = yMin + pi/2; end
    ylim([yMin yMax]);

    yticks(yMin:pi/2:yMax);


    %% ------------------------------------------------------------
    % Generate y-axis labels automatically
    %% ------------------------------------------------------------

    yTickValues = yMin:pi/2:yMax;

    yTickLabels = cell(size(yTickValues));

    for ii = 1:length(yTickValues)

        n = round(2*yTickValues(ii)/pi);

        if n == 0

            yTickLabels{ii} = '0';

        elseif n == 1

            yTickLabels{ii} = '\pi/2';

        elseif n == -1

            yTickLabels{ii} = '-\pi/2';

        elseif n == 2

            yTickLabels{ii} = '\pi';

        elseif n == -2

            yTickLabels{ii} = '-\pi';

        elseif n == 3

            yTickLabels{ii} = '3\pi/2';

        elseif n == -3

            yTickLabels{ii} = '-3\pi/2';

        elseif n == 4

            yTickLabels{ii} = '2\pi';

        elseif n == -4

            yTickLabels{ii} = '-2\pi';

        else

            yTickLabels{ii} = ...
                sprintf('%d\\pi/2',n);

        end

    end

    yticklabels(yTickLabels);

    ylabel('\phi(\theta)');


    %% ------------------------------------------------------------
    % Title
    %% ------------------------------------------------------------

    title(sprintf( ...
        'Polarization winding around BIC, q = %.3f', ...
        qBIC));

end