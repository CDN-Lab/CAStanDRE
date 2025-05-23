function [] = jsonwrite(var,filename)
%used to write a structured variable to a json file

%turn the data matrix into its json equivalent
jsf = jsonencode(var);

%prepare the file name
filesave = strcat(filename,'.json');

%begin a file with writing permission
fid = fopen(filesave,'w');

%write to text file
fprintf(fid,'%s',jsf);

%close the file
fclose(fid);

end