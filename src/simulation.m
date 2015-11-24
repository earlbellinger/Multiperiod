%% Simulation
num_points = 20;
err = 0.1;
x_max = 10*2*pi;
t = 0:0.1:x_max;
p2 = 0.71;
m = sin(t) + 0.4 * sin(2*t+pi/4) ...
    + 0.2 * sin(p2 * t + pi) ...
    + 0.1 * sin(p2 * 8 * t + pi/15);

ust = sort(unifrnd(0, x_max, num_points, 1));
noise = unifrnd(0, err, num_points, 1);
usm = sin(ust) + 0.4*sin(2*ust+pi/4) ...
      + 0.2 * sin(p2 * ust + pi) ...
      + 0.1 * sin(p2 * 8 * ust+pi/15) ...
      + noise;

num_p1 = 4;
num_p2 = 16;
X = ones(length(ust), 1+num_p1+num_p2);
for kk = 2:2:num_p1
    X(:,kk)   = sin(kk/2*ust);
    X(:,kk+1) = cos(kk/2*ust);
end
for kk = 2:num_p2
    X(:,kk+num_p1)   = sin(p2 * kk/2*ust);
    X(:,kk+num_p1+1) = cos(p2 * kk/2*ust);
end

w = (X' * X) \ X' * usm;

wm = ones(1, length(t)) * w(1);
for kk = 2:2:num_p1
    wm = wm + w(kk)   * sin(kk/2*t) ...
            + w(kk+1) * cos(kk/2*t);
end
for kk = 2:2:num_p2
    wm = wm + w(kk+num_p1)   * sin(p2 * kk/2*t) ...
            + w(kk+num_p1+1) * cos(p2 * kk/2*t);
end

[w, FitInfo] = lasso(X, usm);%, 'cv', 10);
lm = ones(1, length(t)) * FitInfo.Intercept(1);
for kk = 2:2:num_p1
    lm = lm + w(kk)   * sin(kk/2*t) ...;
            + w(kk+1) * cos(kk/2*t);
end
for kk = 2:2:num_p2
    lm = lm + w(kk+num_p1)   * sin(p2*kk/2*t) ...
            + w(kk+num_p1+1) * cos(p2*kk/2*t);
end

ymax = max([max(m), max(wm), max(usm), max(lm)]) + 4*err;
ymin = min([min(m), min(wm), min(usm), min(lm)]) - 4*err;

% Plot time series
figure;
plot(t, m, 'color', [.5 0 0], 'LineWidth', 1.5);
xlabel('t');
ylabel('m(t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca, 'box', 'off');
set(gca, 'ylim', [ymin, ymax]);
set(gca, 'xlim', [0, x_max]);
matlab2tikz('mpo.tikz', 'height', '\figureheight', ...
                         'width', '\figurewidth');

% Plot points along curve
figure;
plot(t, m, '--', 'color', [0.5 0.5 0.5], 'LineWidth', 1.5);
xlabel('t');
ylabel('m(t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca, 'box', 'off');
set(gca, 'ylim', [ymin, ymax]);
set(gca, 'xlim', [0, x_max]);
hold on;
errorbar(ust, usm, ones(1,num_points)*2*err, '.', 'MarkerSize', 7.5,...
         'color', [.5 0 0]);
hold off;
matlab2tikz('mpo-points.tikz', 'height', '\figureheight', ...
                               'width',  '\figurewidth');

% Plot OLS
figure;
plot(t, m, '--', 'color', [0.5 0.5 0.5], 'LineWidth', 1);
hold on;
plot(t, wm, '-', 'color', [0 0 0], 'LineWidth', 0.5);
xlabel('t');
ylabel('m(t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca, 'box', 'off');
set(gca, 'ylim', [ymin, ymax]);
set(gca, 'xlim', [0, x_max]);
errorbar(ust, usm, ones(1,num_points)*2*err, '.', 'MarkerSize', 7.5,...
         'color', [.5 0 0]);
hold off;
matlab2tikz('mpo-badfit.tikz', 'height', '\figureheight', ...
                               'width',  '\figurewidth');

% Plot LASSO
figure;
plot(t, m, '--', 'color', [0.5 0.5 0.5], 'LineWidth', 1);
hold on;
plot(t, lm, '-', 'color', [0 0 0], 'LineWidth', 0.5);
xlabel('t');
ylabel('m(t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca, 'box', 'off');
set(gca, 'ylim', [ymin, ymax]);
set(gca, 'xlim', [0, x_max]);
errorbar(ust, usm, ones(1,num_points)*2*err, '.', 'MarkerSize', 7.5,...
         'color', [.5 0 0]);
hold off;
matlab2tikz('mpo-lasso.tikz', 'height', '\figureheight', ...
                              'width', '\figurewidth');
