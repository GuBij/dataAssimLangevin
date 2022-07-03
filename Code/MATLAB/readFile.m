function [data,header]=readFile(file,headerBool,nocol)

 format='';
 for i=1:nocol
   format=strcat(format,'%f');
 end

 fid=fopen(file,'r');
 if headerBool
    header=textscan(fid,'%s[^\n]','Delimiter',{',','"','\t','\b'},'CommentStyle',{'/*','*/'},'MultipleDelimsAsOne',1);
 else
    header='';
 end
 fileInCell=textscan(fid,format,'Delimiter',{'\b','\t',' '},'CommentStyle','Exp','MultipleDelimsAsOne',1);
 fclose(fid);

 data=zeros(length(fileInCell{1}),nocol);
 for i=1:nocol
   data(:,i)=fileInCell{i};
 end

end
