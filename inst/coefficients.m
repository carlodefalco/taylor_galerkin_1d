testnum       = 2;
A             = 1;
omega         = 2;
N             = 200;
L             = 2;
T             = 1;
safety_factor = 1.01;
M             = T/L * N * abs (A) * safety_factor * sqrt (3);
h             = L/N;
dt            = T/M;


