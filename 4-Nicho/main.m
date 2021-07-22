Leroy = importfile("E:\Studium\M2-SeafloorGeodesy\Master2-SeafloorGeodesyUebung\4-Nicho\2301-2303-Leroy.dat", [2, Inf]);

p = Leroy{:,10};
figure
plot(p);

% from data I assume T = 1 year, the data frequency is 2 hour,
% which makes T = 365 * 24 / 2 
T = 365 * 24 / 2;
T = 365 * 24;

% P = a * t + b + c * cos(t/T * 2*pi) + s * sin(t/T * 2 * pi)
t = 1:length(p);
t = t';
A = [t, ones(length(p),1), cos(t/T * 2*pi), sin(t/T * 2 * pi)];

% Ausgleichung
x = (A' * A) \ A' * p;

% The Result realy depends on how you choose T
p_new = A * x;
hold on 
plot(p_new)