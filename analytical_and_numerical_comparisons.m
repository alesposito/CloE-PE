%V4
close all
clear all

% consider hexagonal pillars
s0 = 6;

% check the for kn neighbours
kn = 100;

% check the first ln power decays
ln = 10;

% check the first kc values
kcv = 10*1.6681.^-(0:ln-1);

% analytical expressions
o_3d = @(l)6*zeta(l-2)+2*zeta(l);
o_tl = @(l)12*zeta(l-1)-8*zeta(l)+4/3;
o_db = @(kc,l)exp(1./kc).*(2*polylog(l,exp(-1./kc)) + 6*polylog(l-2,exp(-1./kc)) );
o_3de = @(kc)exp(1./kc).* (8*exp(2./kc)+10*exp(1./kc)+2) ./ (exp(1./kc)-1).^3;

%% testing zeta function approximations for power function

% init
omega3D = []; % Omega for a 3D tissue (power function)
omegaTL = []; % Omega for a 3-layer tissue (power function)
omega3De = []; % Omega for a 3D tissue (exponential decay)
omegaDB = []; % Omega for doubly decaying gradients (3D)


k = (1:kn);

for l=1:ln % power function 
   omega3D(l)  = sum(      (s0*k.^2 +2)./ (k.^l));
   omega3De(l) = sum(      (s0*k.^2 +2) .* exp(-(k-1)./kcv(l))  );  
   omegaTL(l)  = (2/3) *(2+sum(s0*(3*k-2)./ (k.^l)));
   
   for kc=1:ln
       omegaDB(kc,l)  = sum(      (s0*k.^2 +2) .* exp(-(k-1)./kcv(kc)) ./ (k.^l));
   end
end

l = (1:ln);

%% display
figure
subplot(1,3,1)
plot(l,o_3d(l),'-r')
hold on
plot(l,o_tl(l),'-b')
plot(l,omega3D,'or')
plot(l,omegaTL,'ob')
legend({'3d an','tl an','3d num','tl num'})
xlabel('power decay (l)')
ylabel('Omega')
set(gca,'yscale','log','xscale','log','ylim',[1 1e5],'xlim',[0.9 11])
axis square

subplot(1,3,2)
plot(1./kcv,o_3de(kcv),'-r');
hold on
plot(1./kcv,omega3De,'or')
legend({'3de an','3de num'})
xlabel('decay constant (kc^-1)')
ylabel('Omega')
set(gca,'yscale','log','xscale','log','ylim',[1 1e5],'xlim',[0.09 11])
axis square

subplot(1,3,3)
hold on
for kc=1:ln
    plot(l,o_db(kcv(kc),l),'color',[kc./ln 0 0]);
    plot(l,omegaDB(kc,l),'.b');
    str{2*kc-1} = [num2str(kc)  'an'];
    str{2*kc} = [num2str(kc)  'num'];
end

legend(str)
xlabel('power decay (l)')
ylabel('Omega')
set(gca,'yscale','log','xscale','log','ylim',[1 1e3],'xlim',[0.9 11])
axis square




