pattern = [zeros(1000,1)' 0.6*ones(1000,1)'];
period = length(pattern);
I0 = repmat(pattern,20);
I0 = I0(1,:);
p0=zeros(size(I0));
pixsize = (7/period)*1e-06;
lambda = lambda_from_E(20);
z = (0.01:0.01:4)*(7e-06)^2/(4*lambda);

[I_z,~] = fresnel_propagator(I0, p0, pixsize, z, lambda);
carpet = [I0' I_z'];
temp = round(size(carpet,1)/2);
figure, imagesc(carpet(temp:(temp+3*period),:)), colormap gray


pattern = [zeros(1000,1)' pi*ones(1000,1)'];
period = length(pattern);
p0 = repmat(pattern,20);
p0 = p0(1,:);
I0=ones(size(p0));
pixsize = (7/period)*1e-06;
lambda = lambda_from_E(20);
z = (0.01:0.01:4)*(7e-06)^2/(4*lambda);

[I_z_phase,~] = fresnel_propagator(I0, p0, pixsize, z, lambda);
carpet_phase = [I0' I_z_phase'];
temp = round(size(carpet_phase,1)/2);
figure, imagesc(carpet_phase(temp:(temp+3*period),:)), colormap gray



pattern = [zeros(1000,1)' pi*ones(1000,1)'];
period = length(pattern);
p0 = repmat(pattern,30);
p0 = p0(1,:);
pattern2 = [zeros(1000,1)' 0.4*ones(1000,1)'];
I0 = repmat(pattern2,30);
I0=I0(1,:);
pixsize = (7/period)*1e-06;
lambda = lambda_from_E(20);
z = (0.01:0.01:8)*(7e-06)^2/(4*lambda);

[I_z_phase,~] = fresnel_propagator(I0, p0, pixsize, z, lambda);
carpet_phase2 = [I0' I_z_phase'];
temp = round(size(carpet_phase2,1)/2);
figure, imagesc(carpet_phase2(temp:(temp+3*period),:)), colormap gray
