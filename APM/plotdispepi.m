% generates figure 4 of the paper
% plots dispersion relations with epsilon eigenvalue

cyl = CylinderGeometry;

% background properties
cyl.ep = 1; cyl.mu = 1;

% field properties
cyl.k = 1; cyl.a = 1;

% cylinder properties
cyl.mui = 1;

% evaluation properties
morders = 0:2;
norders = 3;
betas = linspace(0.0, 2.0, 500);
betas(betas == cyl.k*cyl.a) = [];

% calculate
epis = zeros([length(morders) length(betas) norders]);
sings = zeros([length(morders) length(betas) 4]);
for k = 1:length(betas)
  cyl.beta = betas(k);
  for l = 1:length(morders)
    cyl.orders = morders(l);
    [roots, poles] = disprootsepi(cyl, norders);
    epis(l,k,:) = permute(roots, [3 1 2]);
    sings(l,k,:) = permute(poles, [3 1 2]);
  end
end

% plot properties
plotm = 0:2;
plotn = 1:3;
radm = betas < cyl.k*cyl.a; 
boundm = betas > cyl.k*cyl.a;

for k = 1:length(plotm)
figure(k)
clf
hold on

% plot roots
for l = 1:length(plotn)
cline(betas(radm), real(epis(plotm(k)+1,radm,plotn(l))), imag(epis(plotm(k)+1,radm,plotn(l))))
cline(betas(boundm), real(epis(plotm(k)+1,boundm,plotn(l))), imag(epis(plotm(k)+1,boundm,plotn(l))))
end 

% light line
ll = line([1 1], [-20 20], 'Color', [0.7 0.7 0.7]);
set(ll, 'LineWidth', 1)

% plot singularities
plot(betas, sings(plotm(k)+1,:,1), 'Color', [0.7 0.7 0.7])
if k ~= 1
  plot(betas, betas.^2./cyl.k^2, 'Color', [0.7 0.7 0.7])
end

%ylim(intellylim(real(epis(plotm(k)+1,:,plotn(l)))))
ylim([-20 20])
caxis([-2 0])
box

title(['m = ' num2str(plotm(k))])
ylabel('\epsilon_c')
xlabel('\beta')

colorbar
end
