% generates figure 3 of the paper
% plots dispersion relations with beta squared eigenvalue

cyl = CylinderGeometry;

% background properties
cyl.ep = 1; cyl.mu = 1;

% field properties
cyl.k = 1; cyl.a = 1;

% cylinder properties
cyl.epi = 12+1i;cyl.mui = 1;

% evaluation properties
morders = 0:4;
norders = 4;
ks = linspace(0.5, 2.0, 100);

% calculate
betas = zeros([length(morders) length(ks) norders]);
sings = zeros([length(morders) length(ks) 4]);
for k = 1:length(ks)
  cyl.k = ks(k);
  for l = 1:length(morders)
    cyl.orders = morders(l);
    [roots, poles] = disprootsbeta(cyl, norders);
    betas(l,k,:) = permute(roots, [3 1 2]);
    sings(l,k,:) = permute(poles(1:4), [3 1 2]);
  end
end

%% plot properties
plotm = morders;
plotn = 1:3;

% plot roots
for k = 1:length(plotn)
figure(k)
clf
hold on

for l = 1:length(plotm)
% boundm = real(betas(l,:,k)) > cyl.ep*cyl.mu.*ks.^2;
  boundm = real(sqrt(betas(l,:,k))) > sqrt(cyl.ep*cyl.mu).*ks & imag(betas(l,:,k)) > 0;
  cline(real(sqrt(betas(l,boundm,k))), ks(boundm), imag(sqrt(betas(1,boundm,k))))
end

% light line
ll1 = line(sqrt(cyl.ep*cyl.mu).*[min(ks) max(ks)], [min(ks) max(ks)], 'Color', [0.7 0.7 0.7]);
ll2 = line(real(sqrt(cyl.epi*cyl.mui)).*[min(ks) max(ks)], [min(ks) max(ks)], 'Color', [0.7 0.7 0.7]);
set(ll1, 'LineWidth', 1)
set(ll2, 'LineWidth', 1)

caxis([0 0.4])
box

title(['n = ' num2str(plotm(k)+1)])
ylabel('k')
xlabel('Re(\beta)')

colorbar
end 
