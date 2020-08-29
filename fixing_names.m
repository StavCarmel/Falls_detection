function fixed=fixing_names(path_of_folder,new_path)
fixed=0;
cd =new_path;
fileList = dir(new_path);
numberOfFiles = length(fileList);


% correcting the file names
for file=1:numberOfFiles
    oldname=fileList(file).name;
    newname=oldname;
    newname(regexp(newname,'\s'))=[];
    newname(regexp(newname,'__'))=[];
    newname((regexp(newname,'v.csv'))+1:(regexp(newname,'v.csv'))+4)=[];
    newname(regexp(newname,'gyp')+2)='r';
    if ~isempty(regexp(newname,'_\d_'))
        newname=insertAfter(newname,regexp(newname,'_\d_'),'0');
    end
    newname(regexp(newname,'cvs')+1)='s';
    newname(regexp(newname,'css')+2)='v';
    if ~isempty(regexp(newname,'3F2F_\d_'))
        newname=insertAfter(newname,regexp(newname,'3F2F_\d_'),'1');
    end
    if ~isempty(strfind(newname,'3BE6'))
        if ~isempty(strfind(newname,'('))
           newname=strcat(newname(1:end-3),'txt'); % so it wont take it for csv
        end
    end
    if ~isempty(strfind(newname,'772E'))
        if ~isempty(strfind(newname,'('))
           newname=strcat(newname(1:end-3),'txt'); % so it wont take it for csv
        end
    end
    if length(newname)<22
        newname(regexp(newname,'67ED_29_08_06_bar')+14)='g';
        newname(regexp(newname,'67ED_29_08_06_gar')+15)='y';
    end
    newname(regexp(newname,'FA5D')+2)='D';
    newname(regexp(newname,'FADD')+3)='5';
    if ~strcmp(oldname,newname)
      f = fullfile(path_of_folder, newname);
      g = fullfile(path_of_folder, oldname);
      movefile(g,f);        %# Rename the fil
    end
end

fixed=1;
end