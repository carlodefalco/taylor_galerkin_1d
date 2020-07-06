testnum       = 5;
A             = 1;
omega         = 2;
N             = 200; % space
L             = 2;
T             = 1;
safety_factor = 1.01;
 M             = T/L * N * abs (A) * safety_factor * sqrt (3); % time
h             = L/N;
dt            = T/(M);

% only N>M
