function [out_oneBig,out_allBig,para]= outlier_detector(delta_xyz)

margin = 0.05; % fraction of extreme data ignored to fit the normal
% distribution to the diplacements.

dist_diff_cut = 0.5; % cutoff of the log difference between the normal
% distribution and the real one.

natoms = size(delta_xyz,1);


out_oneBig = zeros(natoms,1);
out_allBig = zeros(natoms,1);

para = zeros(1,12); %xim, xmax, xmu, xsigma....


idx_outliers = cell(3,1);

for jj=1:3
    
    % 1. histogram of the data
    data = delta_xyz(:,jj);
    [y,x] = hist(data,1e3);
    ynorm = y/trapz(x,y);
    
    % 2. fit a distribution form the margin to 1-margin percentiles
    % 2.1. trim data
    [prov,~] = sort(data);
    data_trimmed =  prov(round(margin*size(data,1)):round((1-margin)*size(data,1)));
    % 2.2. fit a normal distribution to the trimmed data
    PD = fitdist(data_trimmed,'Normal');
    yfit = normpdf(x,PD.mu,PD.sigma);
    
    % 3. calculate the outliers
    % 3.1. calculate the log-distance between the actual data and
    %      the fitted distribution
    res = abs(log(smooth(ynorm,10)')-log(yfit));
    [~,idx] = find(res<dist_diff_cut);
    data_min = min(x(idx));
    data_max = max(x(idx));
    % 3.2. find the indexes of the atoms that are outliers
    [prov,~] = find(data<data_min | data>data_max);
    idx_outliers{jj} = prov;
    % 3.3. save the fitting parameters
    para(1,4*(jj-1)+1) = data_min;
    para(1,4*(jj-1)+2) = data_max;
    para(1,4*(jj-1)+3) = PD.mu;
    para(1,4*(jj-1)+4) = PD.sigma;
    
end

% 4. Select the atoms with at least one large displacement component
%    (out_oneBig) or with all large components (out_allBig)


dxo = zeros(8000,1);
dxo(idx_outliers{1,1})=1;

dyo = zeros(8000,1);
dyo(idx_outliers{2,1})=1;

dzo = zeros(8000,1);
dzo(idx_outliers{3,1})=1;

dxyz = [dxo dyo dzo];

out_oneBig = zeros(8000,1);
out_allBig = zeros(8000,1);
A = sum(dxyz,2);
[idx2,~] = find(A>0);
out_oneBig(idx2) = 1;
[idx3,~] = find(A==3);
out_allBig(idx3) = 1;

% all_outliers = [idx_outliers{1}; idx_outliers{2}; idx_outliers{3}];
% [unique_outliers, w] = unique( all_outliers, 'stable' );
% repeated_outliers = all_outliers(setdiff(1:numel(all_outliers),w));
%
% out_oneBig(setdiff(unique_outliers,repeated_outliers)) = 1;
% out_allBig(repeated_outliers) = 1;


end
