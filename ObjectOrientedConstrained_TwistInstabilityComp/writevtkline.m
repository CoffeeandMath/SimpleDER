function writevtkline(filename,X,Y,Z,varargin)

prec = 7;
fid = fopen(filename, 'w');
%% Generating the Header
fprintf(fid, '# vtk DataFile Version 2.0\n');
%fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'Written badly by Kevin Korner\n');
fprintf(fid, 'ASCII\n');

%% Writing the dataset info (X,Y,Z)
N =  max(size(X));
fprintf(fid, 'DATASET POLYDATA\n');
fprintf(fid, ['POINTS ' num2str(N) ' FLOAT\n']);

for i = 1:N
    fprintf(fid, [num2str(X(i),prec) ' ' num2str(Y(i),prec) ' ' num2str(Z(i),prec) '\n']);
end

fprintf(fid, '\n');

%% Writing the connectivity

fprintf(fid, ['LINES ' num2str(N-1) ' ' num2str(3*(N-1)) '\n']);

for i = 1:(N-1)
    fprintf(fid, [num2str(2) ' ' num2str(i-1) ' ' num2str(i) '\n']);
end


%% Adding other parameters (only scalars at this point)
Natr = max(size(varargin))/2;
if Natr > 0
    
    for i = 1:Natr
        fprintf(fid,'\n\n');
        atr_name = varargin{2*i-1};
        data = varargin{2*i};
        Ndata = max(size(data));
        if i == 1
            fprintf(fid, ['POINT_DATA ' num2str(Ndata) '\n']);
        end
        if ~iscell(data)
            fprintf(fid, ['SCALARS ' atr_name ' FLOAT 1\n']);
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for j = 1:Ndata
                fprintf(fid, [num2str(data(j),prec) '\n']);
            end
        else
            fprintf(fid, ['VECTORS ' atr_name ' FLOAT\n']);
            %fprintf(fid, 'LOOKUP_TABLE default\n');
            for j = 1:Ndata
                fprintf(fid, [num2str(data{j}(1),prec) ' ' num2str(data{j}(2),prec) ' ' num2str(data{j}(3),prec) '\n']);
            end
            
        end
    end
end

%% Ending and Saving File
fclose(fid);
end

