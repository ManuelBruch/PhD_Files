function comma2point_overwrite(FileLocation,FileNameCsv)
%comma2point_overwrite replaces all occurences of comma (",") with point (".") in a text-file.
%   Note that the file is overwritten, which is the price for high speed.

% file    =   memmapfile([FileLocation FileNameCsv],'writable',true);
% dollar  =   uint8('$');
% comma   =   uint8(',');
% point   =   uint8('.');
% file.Data(transpose(file.Data==point))  =   dollar;
% file.Data(transpose(file.Data==comma))  =   point;
% file.Data(transpose(file.Data==dollar)) =   '';

Data = fileread([FileLocation FileNameCsv]);
Data = strrep(Data, '.', '$');
Data = strrep(Data, ',', '.');
Data = strrep(Data, '$', '');
FID = fopen([FileLocation FileNameCsv], 'w');
fwrite(FID, Data, 'char');
fclose(FID);

end

