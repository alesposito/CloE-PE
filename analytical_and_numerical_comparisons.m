close all
clear all

% consider hexagonal pillars
s0 = 6;

% check the for kn neighbours
kn = 100;

% check the first ln power decays
ln = 10;


%% testing zeta function approximations for power function

% init
omega3D = []; % Omega for a 3D tissue
omegaTL = []; % Omega for a 3-layer tissue

k = (1:kn);

for l=1:ln
   omega3D(l) = sum(      (s0*k.^2 +2)./ (k.^l));
   omegaTL(l) = (2/3) *(2+sum(s0*(3*k-2)./ (k.^l)));
end

l = (1:ln);

%
figure
plot(6*zeta(l-2)+2*zeta(l),'-r')
hold on
plot(12*zeta(l-1)-8*zeta(l)+4/3,'-b')
plot(omega3D,'or')
plot(omegaTL,'ob')
legend({'3d zeta','tl zeta','3d num','tl num'})
xlabel('power decay (l)')
ylabel('Omega')
set(gca,'yscale','log','xscale','lin','ylim',[1 1000],'xlim',[0 11])
axis square

%% testing zeta function approximations for exponential function

% init
omega3De = []; % Omega for a 3D tissue

k = (1:kn);

for l=1:ln
   omega3De(l) = sum(      (s0*k.^2 +2) .* exp(-(k-1).*l)  );
end

l = (1:ln);

%
figure
plot(exp(l).* (8*exp(2*l)+10*exp(l)+2) ./ (exp(l)-1).^3,'-r')
hold  on
plot(omega3De,'or')
legend({'3d analytical','3D numerical'})
xlabel('power decay (l)')
ylabel('Omega')
set(gca,'yscale','log','xscale','lin','ylim',[1 100],'xlim',[0 11])
axis square
%%
figure

% linear scale, plot pn/pn0 for power function

k=(1:10)
hold all
for l=1:1:4; plot(k.^-l,'r'); end
for l=1:1:4; plot(exp(-(k-1)*l),'k'); end
xlabel('Cell-to-cell distance (k)')
ylabel('Normalized oncogenic gradient (p_n/p_n_0)')
legend({'power-1','power-2','power','power','exp'})
box on
axis square
set(gca,'yscale','log','xscale','lin')
%set(gca,'yscale','log','xscale','log')


