%This file contains a walkthrough for explaining Gaussian Processes to
%any scientist, including statisticians and non statisticians.
%There is a .pdf file that completes this tutorial.

%%  Demo #1 Univariate Normal %%
% in order to generate a univariate random number, use any real numbers \mu
% and \sigma to define the distribution.
mu = 6.2; %the mean
sigma_squared = 2; %sigma^2 is the variance

%randn returns a random number from the standard normal distribution, which
%is a univariate normal with \mu = 0, \sigma = 1
standard_normal_random_number = randn;

%this is converted to an arbitrary univariate gaussian normal by:
normal_random_number = sqrt(sigma_squared) * standard_normal_random_number + mu;

%now, we can plot what's going on:
for i=1:100000
    standard_normal_random_number(i) = randn;
    normal_random_number(i) = sqrt(sigma_squared) * standard_normal_random_number(i) + mu;
end

hist([standard_normal_random_number' normal_random_number'],50)
legend('standard', 'custom')

%% Demo #2 Multivariate Normal %%
%The same trick works for the multivariate normal. \Mu has to be a vector,
%and \Sigma a symetric semidefinite matrix
Mu = [3 ;0]; %the mean

%a matrix is symetric iff M(a,b) = M(b,a), a positive semidefinite iff
% x*M*x'&gt;=0 for any vector x. These properties are satisfied by all matrices
% which are taken from the matlab function gallery('randcorr',n)
Sigma =  [ 1.0000   -0.9195;  -0.9195    1.0000];

% to find A such that A'*A = Sigma, we find the eigendecompositon of Sigma:
% V'*D*V = Sigma
[V,D]=eig(Sigma);
A=V*(D.^(1/2));

%The standard multivariate normal is just a series of independently drawn
%univariate normal random numbers
standard_random_vector = randn(2,1);

%this is converted to an arbitrary multivariate normal by:
normal_random_vector = A * standard_random_vector + Mu;

%now, we can plot what's going on:
for i=1:1000
    standard_random_vector(:,i) = randn(2,1);
    normal_random_vector(:,i) = A * standard_random_vector(:,i) + Mu;
end
plot(standard_random_vector(1,:),standard_random_vector(2,:),'.')
hold on
plot(normal_random_vector(1,:),normal_random_vector(2,:),'r.')
legend('standard', 'custom')

%% A Kernel Matrix demonstration Brownian motion
n=40;
%The brownian motion kernel is defined through its very simple inverse
inverse=2*eye(n);
inverse=inverse-(triu(ones(n),1)-triu(ones(n),2))';
inverse=inverse-(triu(ones(n),1)-triu(ones(n),2));

subplot(1,3,1);imagesc(inverse);title('Brownian motion K^{-1}')
subplot(1,3,2);imagesc(inv(inverse));title('Brownian motion K')

%draw some samples
r=randn(n,10);
subplot(1,3,3);plot(inverse\r);title('10 samples')

%% Demo #3 Gaussian Process %%
% one way of looking at a Gaussian Process is as a Multivariate Normal with
% an infinite number of dimensions. However, in order to model
% relationships between points, we construct the covariance matrix with a
% function that defines the value of the matrix for any pair of real numbers:
sigma_f = 1.1251; %parameter of the squared exponential kernel
l =  0.90441; %parameter of the squared exponential kernel
kernel_function = @(x,x2) sigma_f^2*exp((x-x2)^2/(-2*l^2));

%This is one of many popular kernel functions, the squared exponential
%kernel. It favors smooth functions. (Here, it is defined here as an anonymous
%function handle)

% we can also define an error function, which models the observation noise
sigma_n = 0.1; %known noise on observed data
error_function = @(x,x2) sigma_n^2*(x==x2);
%this is just iid gaussian noise with mean 0 and variance sigma_n^2s

%kernel functions can be added together. Here, we add the error kernel to
%the squared exponential kernel)
k = @(x,x2) kernel_function(x,x2)+error_function(x,x2);

%We can find the mean and the variance of the GP at each point
prediction_x=-2:0.01:1;
for i=1:length(prediction_x)
    mean(i) = 0;
    variance(i) = k(prediction_x(i),prediction_x(i));
end
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,fliplr(lower)],color),'EdgeColor',color);
plot_variance(prediction_x,mean-1.96*variance,mean+1.96*variance,[0.9 0.9 0.9])
hold on
set(plot(prediction_x,mean,'k'),'LineWidth',2)

%% Demo #4 Gaussian Process Sampling %%
% now, we would like to sample from a Gaussian Process defined by this
% kernel. First, for the subset of points we are interested to plot, we
% construct the kernel matrix (using our kernel function)
K=zeros(length(prediction_x),length(prediction_x));
for i=1:length(prediction_x)
    for j=i:length(prediction_x)%We only calculate the top half of the matrix. This is an unnecessary speedup trick
        K(i,j)=k(prediction_x(i),prediction_x(j));
    end
end
K=K+triu(K,1)'; % We can use the upper half of the matrix and copy it to the
%bottom, because it is symetrical

[V,D]=eig(K);
A=V*(D.^(1/2));

%Then, we use the kernel as the covariance matrix of a multivariate normal
clear gaussian_process_sample;
for i=1:7
    standard_random_vector = randn(length(prediction_x),1);
    gaussian_process_sample(:,i) = A * standard_random_vector;
end

plot(prediction_x,real(gaussian_process_sample))

%% Demo #5 Gaussian Process Regression %%
%initialize observations
X_o = [-1.5 -1 -0.75 -0.4 -0.3 0]';
Y_o = [-1.6 -1.3 -0.5 0 0.3 0.6]';

K = zeros(length(X_o));
for i=1:length(X_o)
    for j=1:length(X_o)
        K(i,j)=k(X_o(i),X_o(j));
    end
end

K_ss=zeros(length(prediction_x),length(prediction_x));
for i=1:length(prediction_x)
    for j=i:length(prediction_x)%We only calculate the top half of the matrix. This an unnecessary speedup trick
        K_ss(i,j)=k(prediction_x(i),prediction_x(j));
    end
end
K_ss=K_ss+triu(K_ss,1)'; % We can use the upper half of the matrix and copy it to the

K_s=zeros(length(prediction_x),length(X_o));
for i=1:length(prediction_x)
    for j=1:length(X_o)%We only calculate the top half of the matrix. This an unnecessary speedup trick
        K_s(i,j)=k(prediction_x(i),X_o(j));
    end
end

Mu = (K_s/K)*Y_o;
Sigma = 1.96*sqrt(diag(K_ss-K_s/K*K_s'));

figure
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,fliplr(lower)],color),'EdgeColor',color);
plot_variance(prediction_x,(Mu-Sigma)',(Mu+Sigma)',[0.8 0.8 0.8])
hold on
set(plot(prediction_x,Mu,'k-'),'LineWidth',2)
set(plot(X_o,Y_o,'r.'),'MarkerSize',15)

% this gives a poor model, because we aren't using good parameters to model
% the function. In order to get better parameters, we can maximize evidence
evidence = exp((Y_o'/K*Y_o+log(det(K))+length(Y_o)*log(2*pi))/-2);
title (['this plot has evidence ' num2str(evidence)])

legend('confidence bounds','mean','data points','location','SouthEast')

%% Demo #5.2 Sample from the Gaussian Process posterior
clearvars -except k prediction_x K X_o Y_o

%We can also sample from this posterior, the same way as we sampled before:
K_ss=zeros(length(prediction_x),length(prediction_x));
for i=1:length(prediction_x)
    for j=i:length(prediction_x)%We only calculate the top half of the matrix. This an unnecessary speedup trick
        K_ss(i,j)=k(prediction_x(i),prediction_x(j));
    end
end
K_ss=K_ss+triu(K_ss,1)'; % We can use the upper half of the matrix and copy it to the

K_s=zeros(length(prediction_x),length(X_o));
for i=1:length(prediction_x)
    for j=1:length(X_o)
        K_s(i,j)=k(prediction_x(i),X_o(j));
    end
end

[V,D]=eig(K_ss-K_s/K*K_s');
A=real(V*(D.^(1/2)));

for i=1:7
    standard_random_vector = randn(length(A),1);
    gaussian_process_sample(:,i) = A * standard_random_vector+K_s/K*Y_o;
end
hold on
plot(prediction_x,real(gaussian_process_sample))

set(plot(X_o,Y_o,'r.'),'MarkerSize',20)

%% Demo 6 finding parameters of a Gaussian Process with grid search %%
sigma_n = 0.2;% we know the amount of noise from the data
sigma_range = 0.01:0.04:4;
l_range = 0.01:0.04:2;

evidence=zeros(length(sigma_range),length(l_range));
for i=1:length(sigma_range)
    for j=1:length(l_range)
        evidence(i,j)=evidence_2_param_GP([sigma_range(i) l_range(j)],sigma_n);
    end
end
imagesc(sigma_range,l_range,real(evidence)');colorbar
title('inverse log evidence of the squared exponential kernel')
xlabel('sigma parameter')
ylabel('l parameter')
hold on
%get the max l and sigma
[v location]=min(evidence(:));
l_location=floor(location/length(sigma_range))+1;
sigma_location=mod(location,l_location*length(sigma_range)-length(sigma_range));
l = l_range(l_location);
sigma = sigma_range(sigma_location);

set(plot(sigma,l,'rx'),'MarkerSize',25)

% finding parameters of a Gaussian Process with fminsearch %%
starting_place = randn(2,10);
sigma_n = 0.2;
for i=1:size(starting_place,2)
    x_guess(:,i)=fminsearch(@evidence_2_param_GP,starting_place(:,i),[],sigma_n);
    error(i)=evidence_2_param_GP(x_guess(:,i),sigma_n);
end
[v location]=min(error);

disp(['The optimal evidence is ' num2str(exp(-evidence_2_param_GP(x_guess(:,location),sigma_n))) ','])
disp(['(for parameters \sigma=' num2str(x_guess(1,location)) ' and l=' num2str(x_guess(2,location)) ')'])
x_guess=abs(x_guess);
set(plot(x_guess(1,location),x_guess(2,location),'r.'),'MarkerSize',25)

legend('grid search minimum', 'matlab''s fminsearch minimum')

%% Demo 7 Gaussian Process with input in 2-D
clear
close all
sigma_f=0.5;
l=1;
kernel_function = @(x,x2) sigma_f^2*exp(((x-x2)'*(x-x2))/(-2*l^2));

grid_size=15;
[x1 x2]=meshgrid(1:grid_size,1:grid_size);
prediction_x=[x1(:)';x2(:)']./5;

for i=1:grid_size^2
    for j=1:grid_size^2
        K(i,j)=kernel_function(prediction_x(:,i),prediction_x(:,j));
    end
end

[V,D]=eig(K);
A=V*(D.^(1/2));

for k=1:10
    % prediction_x is a matrix where each collumn is a test location in x
    standard_random_vector = randn(length(prediction_x),1);
    gaussian_process_sample = A * standard_random_vector;
    %plot3(prediction_x(1,:),prediction_x(2,:),real(gaussian_process_sample(:,i)));
    %grid on
    surf(x1,x2,reshape(real(gaussian_process_sample),grid_size,grid_size));
    xlabel('x_1')
    ylabel('x_2')
    zlabel('y')
    pause
end

%% Demo 8 (offline) Flying through a Guassian Process sample in 2-D
clear
close all
sigma_f=0.5;
l=1;
kernel_function = @(x,x2) sigma_f^2*exp(dot(x-x2,x-x2)/(-2*l^2));

grid_size=20;
animation_length=100;
tic
[x1 x2]=meshgrid(1:grid_size,1:grid_size+animation_length);
toc
% prediction_x is a matrix where each collumn is a test location in x
prediction_x=[x1(:)';x2(:)']'./5;
% for i=1:length(prediction_x)
%     for j=1:length(prediction_x)
%         K(i,j)=kernel_function(prediction_x(:,i),prediction_x(:,j));
%     end
% end
%new version
tic
n=size(prediction_x,1);
K=prediction_x*prediction_x'/sigma_f^2;
d=diag(K);
K=K-ones(n,1)*d'/2;
K=K-d*ones(1,n)/2;
K=exp(K);
toc

[V,D]=eig(K);
A=V*(D.^(1/2));

standard_random_vector = randn(length(prediction_x),1);
gaussian_process_sample = real(A * standard_random_vector)./3;

sample_rect=reshape(real(gaussian_process_sample),grid_size+animation_length,grid_size);
i=1;
surf(x1(i:i+grid_size-1,:),x2(i:i+grid_size-1,:),sample_rect(i:i+grid_size-1,:));

for i=1:animation_length
    %plot3(prediction_x(1,:),prediction_x(2,:),real(gaussian_process_sample(:,i)));
    %grid on
    surf(x1(i:i+grid_size-1,:),x2(i:i+grid_size-1,:),sample_rect(i:i+grid_size-1,:));
    axis([0 20 i i+grid_size -1 1])
    xlabel('x_1')
    ylabel('x_2')
    zlabel('y')
    pause(0.1)
end

%% Demo 9 (online) Flying through a Guassian Process sample in 2-D
clear
close all
sigma_f=7;

grid_size=50;
[x1 x2]=meshgrid(1:grid_size);
x1_place=grid_size;
% prediction_x is a matrix where each collumn is a test location in x
prediction_x=[x1(:)';x2(:)']';
% for i=1:length(prediction_x)
%     for j=1:length(prediction_x)
%         K(i,j)=kernel_function(prediction_x(:,i),prediction_x(:,j));
%     end
% end
%new version
K=rbf(prediction_x,sigma_f);
[V,D]=eig(K);
A=V*(D.^(1/2));

standard_random_vector = randn(length(prediction_x),1);
gaussian_process_sample = real(A * standard_random_vector);

sample_rect=reshape(real(gaussian_process_sample),grid_size,grid_size);
surf(x1,x2,sample_rect);
i=1;
axis([1 grid_size i i+grid_size-1 -2 2])
xlabel('x_1')
ylabel('x_2')
zlabel('y')

for i=1:2000
    %generate next x coordinates
    x1_place = x1_place + 1;
    [x1_next x2_next]=meshgrid(x1_place,1:grid_size);
    x1_old=x1(:,end-2:end);
    x2_old=x2(:,end-2:end);
    x1=[x1(:,2:end) x1_next];
    x2=[x2(:,2:end) x2_next];
    previous_x=[x1_old(:) x2_old(:)];
    new_x=[x1_next(:) x2_next(:)];

    %generate y for these coordinates, given already sampled y. To reduce
    %complexity, we only look at the y of the previous frame, not all observed y
    tic
    K=rbf(previous_x,sigma_f)+eye(length(previous_x))*10^-5;
    K_ss=rbf(new_x,sigma_f)+eye(length(new_x))*10^-5;
    K_s=rbf(new_x,sigma_f,previous_x);

    [V,D]=eig(K_ss-K_s/K*K_s');
    A=real(V*(D.^(1/2)));

    standard_random_vector = randn(length(K_ss),1);
    %gaussian_process_sample
    %real(gaussian_process_sample(grid_size+1:end))
    new_sample = A * standard_random_vector+K_s/K*real(gaussian_process_sample(end-length(previous_x)+1:end));
    gaussian_process_sample = [gaussian_process_sample(length(new_sample)+1:end);new_sample];
    toc
    sample_rect=reshape(real(gaussian_process_sample),grid_size,grid_size);

    surf(x1,x2,sample_rect);
    axis([i i+grid_size-1 1 grid_size -2 2])
    xlabel('x_1')
    ylabel('x_2')
    zlabel('y')

    pause(0.05)
end

%% Demo 10 Sampling the Gaussian Process prior in 3D
%the previously used kernel_function and error_function can be extended to
%support vector input. This is an natural case of the squared exponential
%kernel, with the same parameters and properties. (note that this is the
%rotationally invariant version)
kernel_function_m = @(x,x2) sigma_f^2*exp((x-x2)'*(x-x2)/(-2*l^2));
%here, the error function needs to test whether two vectors are exactly the
%same. This is done by counting the number of matches between x and x2
error_function_m = @(x,x2) sigma_n^2*(sum(x==x2)==length(x));
%k_m is used instead of k
k_m = @(x,x2) kernel_function_m(x,x2)+error_function_m(x,x2); 

%the resolution option here allows to change the granularity of the sampled
%Gaussian Process. Since resolution^2 samples need to be taken, this cannot
%be too high. In even higher dimensions, taking samples from a Gaussian
%Process at every point of a grid becomes prohibitive, so they cannot be
%visualised by their samples.
resolution=5;
%generate the grid where to take samples, and save it the same way as
%prediciton_x in previous Demos.
[a b]=meshgrid(linspace(0,1,resolution));
%Before, prediciton_x had n*1 dimensions, now it will have n*2
prediction_x=[a(:) b(:)]';

%This is done exactly as in Demo 4
K=zeros(size(prediction_x,2),size(prediction_x,2));
for i=1:size(prediction_x,2)
    for j=i:size(prediction_x,2)%We only calculate the top half of the matrix.
        %(This is an unnecessary speedup trick)

        %here, we use the vectors prediction_x(:,i) for each point
        K(i,j)=k_m(prediction_x(:,i),prediction_x(:,j));
    end
end
K=K+triu(K,1)'; % We can use the upper half of the matrix and copy it to the
%bottom, because it is symetrical

[V,D]=eig(K);
A=V*(D.^(1/2));

%take a single sample, and plot the surface. Note that this becomes
%impossible in higher dimensions.
standard_random_vector = randn(size(prediction_x,2),1);
sample=A * standard_random_vector;
surf(reshape(sample,resolution,resolution))

%an important note: a Gaussian Process with 1D input is a continuous line, but for a
%2D input, it is a continuous surface. That's because the Gaussian Process defines a
%value for every possible point everywhere in the x-y space.

%% Demo 11 Sampling the Gaussian Process posterior in many dimensions
dimensions=4;
%this initialises the random number generator to create the same numbers
%each time it is executed, for repeatability.
rng('default');
%generate 7 random observations in 4 dimensions
X_o = rand(dimensions,7);
Y_o = [-2 -1.6 -1.3 -0.5 0 0.3 0.6]';

K = zeros(size(X_o,2));
for i=1:size(X_o,2)
    for j=1:size(X_o,2)
        K(i,j)=k_m(X_o(:,i),X_o(:,j));
    end
end

prediction_x=[0 0 0 0; 1 0 0 0; 1 1 1 1]';

K_ss=zeros(size(prediction_x,2));
for i=1:size(prediction_x,2)
    for j=i:size(prediction_x,2)%We only calculate the top half of the matrix. This an unnecessary speedup trick
        K_ss(i,j)=k_m(prediction_x(:,i),prediction_x(:,j));
    end
end
K_ss=K_ss+triu(K_ss,1)'; % We can use the upper half of the matrix and copy it to the

K_s=zeros(size(prediction_x,2),size(X_o,2));
for i=1:size(prediction_x,2)
    for j=1:size(X_o,2)%We only calculate the top half of the matrix. This an unnecessary speedup trick
        K_s(i,j)=k_m(prediction_x(:,i),X_o(:,j));
    end
end

%calculate Mu and Sigma according to the equation on page 13
Mu = (K_s/K)*Y_o;
Sigma = K_ss-K_s/K*K_s';