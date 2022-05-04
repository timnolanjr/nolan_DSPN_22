%% Generate data for DSPN Final Project
% Tim Nolan, tnolan@andrew.cmu.edu

%Set seed for reproducibility
rng(123)

Fs = 10;
t=0:1/Fs:7200-1/Fs; % 3600s = 1h

%% 1. Make Nosie Data
iid_noise=nirs.testing.simARNoise([],t,0,.33);% set P=0 for IID noise
sin_noise = nirs.testing.simPhysioNoise(iid_noise,.25,1,2);
ar_noise = nirs.testing.simARNoise([],t,5,.33); % AR(5)

sigma_iid = var(reshape(iid_noise.data,[],1));
sigma_ar = var(reshape(ar_noise.data,[],1));
sigma_sin = var(reshape(sin_noise.data,[],1));
sf = sqrt(sigma_ar/sigma_sin);
sf_iid = sqrt(sigma_ar/sigma_iid);

mu = mean(sin_noise.data,'all');
sin_noise.data = mu + sf * (sin_noise.data - mu);

mu_iid = mean(iid_noise.data,'all');
iid_noise.data = mu_iid + sf_iid * (iid_noise.data - mu_iid);

sinar_noise = nirs.testing.simPhysioNoise(ar_noise,.25,1,2);


%% 2. Set Regression slopes
beta = [.1,.1,.5,.5,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,zeros(1,8)];
beta(2:2:end)=-.5*beta(2:2:end);

%SD Pairs with true relationships
channels=[1	1;1	1;2	1;2	1;2	2;2	2;3	2;3	2;3	3;3	3;4	3;4	3;4	4;4	4;5	4;5	4;5	5;5	5;6	5;6	5;6	6;6	6;7	6;7	6];

%% 3. Simulate Data
stimFunc = nirs.testing.blockedStimDesign(t,30,30,1,1);

[rawdata_iid,truth_iid] = nirs.testing.simData(iid_noise,stimFunc,beta,channels);
[rawdata_ar,truth_ar] = nirs.testing.simData(ar_noise,stimFunc,beta,channels);
[rawdata_sin,truth] = nirs.testing.simData(sin_noise,stimFunc,beta,channels);
[rawdata_sinar,truth_sinar] = nirs.testing.simData(sinar_noise,stimFunc,beta,channels);

%% 4. Preprocess
job = nirs.modules.OpticalDensity();
job = nirs.modules.BeerLambertLaw(job);
job = nirs.modules.Resample(job);
job.Fs = 5;

Fs = 5;
t=0:1/Fs:3600-1/Fs;

hbdata_iid = job.run(rawdata_iid);
hbdata_ar = job.run(rawdata_ar);
hbdata_sin = job.run(rawdata_sin);
hbdata_sinar = job.run(rawdata_sinar);

%% 5. Get X data
X=hbdata_ar.getStimMatrix;
basis = nirs.design.basis.Canonical;
h = basis.getFilter(hbdata_ar.Fs);
X0 = filter(h,1,X);


%% 6. Run AR-IRLS to get residuals
job2 = nirs.modules.AR_IRLS();
resid_iid = job2.run(hbdata_iid);
resid_ar = job2.run(hbdata_ar);
resid_sin = job2.run(hbdata_sin);
resid_sinar = job2.run(hbdata_sinar);


%% 7. Save
df = [hbdata_iid.data,X0];
save('./data/y_iid_raw.mat','df','-v6')

df = [resid_iid.ydata,resid_iid.xdata];
save('./data/y_iid_resid.mat','df','-v6')

df = [hbdata_ar.data,X0];
save('./data/y_ar_raw.mat','df','-v6')

df = [resid_ar.ydata,resid_ar.xdata];
save('./data/y_ar_resid.mat','df','-v6')

df = [hbdata_sin.data,X0];
save('./data/y_sin_raw.mat','df','-v6')

df = [resid_sin.ydata,resid_sin.xdata];
save('./data/y_sin_resid.mat','df','-v6')

df = [hbdata_sinar.data,X0];
save('./data/y_sinar_raw.mat','df','-v6')

df = [resid_sinar.ydata,resid_sinar.xdata];
save('./data/y_sinar_resid.mat','df','-v6')

%% 8. Plot
corder=cbrewer('qual','Paired',8);
figure
tiledlayout(5,1)
nexttile
plot(hbdata_iid.time(1:600),hbdata_iid.data(1:600,1),'Color',corder(2,:))
hold on
%plot(hbdata_iid.time(1:600),hbdata_iid.data(1:600,2),'Color',corder(1,:))
plot(hbdata_iid.time(1:600),hbdata_iid.data(1:600,5),'Color',corder(4,:))
%plot(hbdata_iid.time(1:600),hbdata_iid.data(1:600,6),'Color',corder(3,:))
plot(hbdata_iid.time(1:600),hbdata_iid.data(1:600,23),'Color',corder(6,:))
%plot(hbdata_iid.time(1:600),hbdata_iid.data(1:600,24),'Color',corder(5,:))
ax=gca;
ax.YLim=[-100, 100];
xlabel('Time(s)','FontSize',14)
ylabel('Intensity (au)','FontSize',14)
title('IID Noise','FontSize',18)
legend('Ch 1, SNR = 0.1','Ch 5, SNR = 1.0','Ch 23, SNR = 10.','FontSize',14)

nexttile
plot(hbdata_ar.time(1:600),hbdata_ar.data(1:600,1),'Color',corder(2,:))
hold on
%plot(hbdata_ar.time(1:600),hbdata_ar.data(1:600,2),'Color',corder(1,:))
plot(hbdata_ar.time(1:600),hbdata_ar.data(1:600,5),'Color',corder(4,:))
%plot(hbdata_ar.time(1:600),hbdata_ar.data(1:600,6),'Color',corder(3,:))
plot(hbdata_ar.time(1:600),hbdata_ar.data(1:600,23),'Color',corder(6,:))
%plot(hbdata_ar.time(1:600),hbdata_ar.data(1:600,24),'Color',corder(5,:))
ax=gca;
xlabel('Time(s)','FontSize',14)
ylabel('Intensity (au)','FontSize',14)
title('AR(5) Noise','FontSize',18)
legend('Ch 1, SNR = 0.1','Ch 5, SNR = 1.0','Ch 23, SNR = 10.','FontSize',14)

nexttile
plot(hbdata_sin.time(1:600),hbdata_sin.data(1:600,1),'Color',corder(2,:))
hold on
%plot(hbdata_sin.time(1:600),hbdata_sin.data(1:600,2),'Color',corder(1,:))
plot(hbdata_sin.time(1:600),hbdata_sin.data(1:600,5),'Color',corder(4,:))
%plot(hbdata_sin.time(1:600),hbdata_sin.data(1:600,6),'Color',corder(3,:))
plot(hbdata_sin.time(1:600),hbdata_sin.data(1:600,23),'Color',corder(6,:))
%plot(hbdata_sin.time(1:600),hbdata_sin.data(1:600,24),'Color',corder(5,:))
ax=gca;
xlabel('Time(s)','FontSize',14)
ylabel('Intensity (au)','FontSize',14)
title('Sinusoidal Noise','FontSize',18)
legend('Ch 1, SNR = 0.1','Ch 5, SNR = 1.0','Ch 23, SNR = 10.','FontSize',14)

nexttile
plot(hbdata_sinar.time(1:600),hbdata_sinar.data(1:600,1),'Color',corder(2,:))
hold on
%plot(hbdata_sinar.time(1:600),hbdata_sinar.data(1:600,2),'Color',corder(1,:))
plot(hbdata_sinar.time(1:600),hbdata_sinar.data(1:600,5),'Color',corder(4,:))
%plot(hbdata_sinar.time(1:600),hbdata_sinar.data(1:600,6),'Color',corder(3,:))
plot(hbdata_sinar.time(1:600),hbdata_sinar.data(1:600,23),'Color',corder(6,:))
%plot(hbdata_sinar.time(1:600),hbdata_sinar.data(1:600,24),'Color',corder(5,:))
ax=gca;
xlabel('Time(s)','FontSize',14)
ylabel('Intensity (au)','FontSize',14)
title('Physiological Noise','FontSize',18)
legend('Ch 1, SNR = 0.1','Ch 5, SNR = 1.0','Ch 23, SNR = 10.','FontSize',14)

nexttile
plot(hbdata_iid.time(1:600),X(1:600),'Color','k')
hold on
plot(hbdata_iid.time(1:600),X0(1:600),'Color','k','LineStyle','--')
xlabel('Time(s)','FontSize',14)
ylabel('Intensity (au)','FontSize',14)
title('Experimental Design Matrix','FontSize',18)
legend('Design, X','X_{HRF}','FontSize',14)
