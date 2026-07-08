% Test harness for seedingRate.m under GNU Octave
% Uses a small reproducible clone-frequency matrix.

% 4 tumors (rows) x 5 clones (cols). Unnormalized counts on purpose
% so we also exercise the internal row-normalization step.
y = [ 50 30 10  5  5;   % tumor 1 (we'll call this the primary)
      40 35 15  7  3;   % tumor 2
      10 20 40 20 10;   % tumor 3
       5 10 25 30 30 ]; % tumor 4

disp('=== Input matrix y ===');
disp(y);

% --- Case A: primary tumor designated (m = 1), Proposition 7 path ---
disp('=== Case A: m = 1 (primary designated) ===');
[kA, xA, DKLA] = seedingRate(y, 1);
disp('Inferred seeding rates k:'); disp(kA);
disp('Source clone freqs x:');     disp(xA);
disp('KL divergence DKL:');        disp(DKLA);

% --- Case B: primary tumor NOT designated (m = 0), Corollary 7.4 path ---
disp('=== Case B: m = 0 (primary not designated) ===');
[kB, xB, DKLB] = seedingRate(y, 0);
disp('Inferred seeding rates k:'); disp(kB);
disp('Source clone freqs x:');     disp(xB);
disp('KL divergence DKL:');        disp(DKLB);
