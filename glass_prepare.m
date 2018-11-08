% prepare the data

Raw_Data = load('glass.data.txt');
Num_Data_Pnts = size(Raw_Data,1);
Dim = size(Raw_Data,2);

data = Raw_Data;


% Center the data
mu = mean(data);
data = data - ones(Num_Data_Pnts,1) * mu;

% Normalize the data, diveded by std
for d = 1:Dim
    data(:,d) = data(:,d)./std(data(:,d));
end


% Pick the columns for parameters
data = data(:,2:Dim-1);
% Pick the column for dependent variable
dvar = Raw_Data(:,Dim);

save data.parameters.norm.dat data -ascii;
save data.dependent_var.norm.dat dvar -ascii;