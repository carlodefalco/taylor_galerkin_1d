coefficients

x = linspace (0, L, N+1);
t = linspace (0, T, M+1);

for n = 1 : M
  t(n)
  [u, ux, ut] = example_solution (t(n), x);
  plot (x, u);
  drawnow
  pause (.1)
end
