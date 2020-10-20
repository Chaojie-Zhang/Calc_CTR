clear all;
clc;
close all;

% parpool('local',3);
    %parallel program  
    
load('data_46.mat');
x1 = data_46{1} * 0.1273; %um
p1 = data_46{2}; %gamma beta pz
x2 = data_46{3} * 0.1273; %um
p2 = data_46{4};
x3 = data_46{5} * 0.1273; %um
p3 = data_46{6};
q = data_46{7};


figure(1);
plot(x1,p1,'.');
xlabel('x1')
ylabel('p1');

%% Input


mm = 1e-1; %mm to cm
cgs_c=3*10^10; %cm/s
cgs_e= 4.8*10^-10; %4.803 204 27 ¡Á 10^-10 Fr
f = 300 * mm; %Distance of measurement
total_charge = 20e-12 ; %pC
q = q/sum(q) * total_charge;
qn = q/1.6e-19;
x2max = 2/2; %cm
n = 200;   %num of grid
dx = 2*x2max/n;
x = -x2max:dx:x2max; % [-x2max, x2max]
y = x;
nw = 1000; %num of frequency bins
w = 1e13:(5e15-1e13)/nw:5e15;  % frequency range from 1e13 Hz to 5e15 Hz

% %read file
% file_id = fopen('electron_40gamma_0.040mrad_more.o', 'r');
% data = zeros(47257, 17);
% count = 1;
% while feof(file_id) == 0
%    temp = fread(file_id, 17, 'double');
%    if feof(file_id) ~= 0
%        break;
%    end
%    data(count, :) = temp;
%    count = count+1;
% end


num = 63960; % number of particles
num_sample = 10000; %sampled num
idx = randperm(num);
idx = idx(1:num_sample);
bmatrix = zeros(n,n,nw); 
rmatrix = zeros(n,n,n);
% pdata = data(:,[13,14,15]); %px, py,pz MeV
% r1data = data(:,[10,11,12]) * mm; %cm
% energydata = data(:,17); %MeV
% tdata = data(:,16) * 1e-9; %s

pdata = zeros(num_sample, 3);
r1data = zeros(num_sample,3);
pdata(:,3) = p1(idx) * 0.511; % gamma beta * m0C^2
pdata(:,1) = p2(idx) * 0.511; 
pdata(:,2) = p3(idx) * 0.511;
r1data(:,3) = 0.0;          %r1 defines the position of particle on the radiator
% r1data(:,1) = 0.0;
% r1data(:,2) = 0.0;
r1data(:,1) = x2(idx) * 1e-4; %um to cm
r1data(:,2) = x3(idx) * 1e-4; %um to cm
energydata = sqrt(p1.^2 + p2.^2 + p3.^2 + 1) .* 0.511; % total Energy
tdata = x1 * 1e-6 ./ (p1 * 0.511 ./energydata * 3e8); %Here change position in z into time arrived.

%% Calculation


n = length(x);
nw = length(w);
r2matrix = zeros(n,n,3); %r2 defines the position of the receptive screen
r2matrix(:,:,1) = repmat(x,n,1);
r2matrix(:,:,2) = repmat(y',1,n);
r2matrix(:,:,3) = f * ones(n,n);

bmatrix(n,n,nw) = 0;

%tdist = sigmat * randn(1,num);
for i = 1:num_sample
    p = pdata(i,:);
    r12 = r1data(i,:);
    energy = energydata(i);
    t = tdata(i);
    if mod(i,10) == 0
        i
    end
    
    absbeta = (1-(0.511/energy)^2)^0.5;
    beta(1,1,:) = absbeta * p/norm(p);
    r1(1,1,:) = r12;
    
    rmatrix = r2matrix - repmat(r1,size(r2matrix,1),size(r2matrix,2));
    normr = sqrt(sum(rmatrix.^2,3));
    rn = rmatrix./repmat(normr,1,1,3);
    bcn = cross(repmat(beta,n,n),rn,3);
    bdn = dot(repmat(beta,n,n),rn,3);
    for l = 1:nw
%      bmatrix(:,:,l) = cgs_e./(cgs_c * normr) .* ... %original form with
%      same charge
%             sqrt(sum((bcn./repmat((1-bdn),1,1,3) + (-bcn)./repmat(1+bdn,1,1,3)).^2,3))...
%             .* exp(1i * w(l) * (normr/cgs_c + t))+bmatrix(:,:,l);
%         
     bmatrix(:,:,l) = qn(i) * cgs_e./(cgs_c * normr) .* ...
    sqrt(sum((bcn./repmat((1-bdn),1,1,3) + (-bcn)./repmat(1+bdn,1,1,3)).^2,3))...
    .* exp(1i * w(l) * (normr/cgs_c + t))+bmatrix(:,:,l);
    end
end
%bmatrix = qn/num_sample * bmatrix * num/num_sample;
bmatrix =  bmatrix * num/num_sample; % fix sample 
dedw = 2* cgs_c / (2*pi) * abs(bmatrix).^2/(4*pi*1) * dx^2 .* f./repmat(normr,1,1,nw);
dedws = squeeze(sum(sum(dedw,1),2)) * 1e-7 * 1e9 * 1e12; %nJ/THz
dedlambdas = dedws'.*w.^2/(2*pi*cgs_c) * 1e-12 * 1e-4; %nJ/um
lambda = 2*pi*cgs_c./w * 1e4;

%% Figures of the results
%load('2fs_20000_76um_1cm_f646_15d_d2_osiris.mat')
figure(1);
[xx,yy] = meshgrid(x,y);
mesh(xx,yy,sum(dedw,3));
view(2);
xlabel('x');
ylabel('y');

figure(2);
plot(w*1e-12/(2*pi),dedws, '-','linewidth',2,'Color',[255,69,0]/255);
xlabel('\omega(THz)');
ylabel('Energy(nJ/THz)');
%axis([0 1000 0 8e-3]);
set (gcf,'Position',[500,250,500,450], 'color','w','PaperpositionMode','auto');
set(gca,'Fontsize',12 ,'linewidth',1.15,'FontWeight','bold','ticklength',[0.015 0.015]);


figure(3);
plot(lambda(lambda>0 & lambda<15),dedlambdas(lambda>0 & lambda<15),'linewidth',2);
xlabel('\lambda(um)');
ylabel('Energy(nJ/um)');
%ylabel('Energy');

set (gcf,'Position',[500,250,500,450], 'color','w','PaperpositionMode','auto');
set(gca,'Fontsize',12 ,'linewidth',1.15,'FontWeight','bold','ticklength',[0.015 0.015]);



totalenergy = sum(dedws(lambda<8&lambda>2)*(w(2)-w(1))*1e-12);


%save('47257_76um_1cm_f646_15d_d2_osiris_more.mat');

fclose('all');




