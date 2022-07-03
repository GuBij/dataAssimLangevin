function [data,header]=readPGFile(file,headerBool,nocol,arcBool)

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
 data = [];
 if ( arcBool )
   blockTitleFormat1 = '%f%f%s%s%f'; 
   blockTitleFormat2 = '%d%d%d%d%f%f%f%f%f%f%f'; 
   arcNr = 0;
 else
   blockTitleFormat1 = 'Run No. %f                      Post Location (degrees)';
   blockTitleFormat2 = '%s%d%d%d%d%d%d';
 end
 expNr = 0;
 while ( true )
  blockID1=textscan(fid,blockTitleFormat1,1,'Delimiter',{'\b','\t',' ',','},'CommentStyle',{'/*','*/'},'MultipleDelimsAsOne',1);
  if ( feof(fid) )
    break;
  end
  blockID2=textscan(fid,blockTitleFormat2,1,'Delimiter',{'\b','\t',' ',','},'MultipleDelimsAsOne',1);
  expNr = blockID1{1}; 
  if ( arcBool )
   arcNr = blockID1{2};
   noRows = blockID2{1};
   block=zeros(noRows,nocol+2);
   block(:,1:2)=[expNr*ones(noRows,1),arcNr*ones(noRows,1)];
   blockInCell=textscan(fid,format,noRows,'Delimiter',{'\b','\t',' '},'MultipleDelimsAsOne',1);
   for i=1:nocol
    block(:,i+2)=blockInCell{i};
   end
   block=block(blockID2{2}:blockID2{3},:);
  else
   noRows = 9;
   block=zeros(noRows,nocol+1);
   block(:,1)=expNr*ones(noRows,1);
   blockInCell=textscan(fid,format,noRows,'Delimiter',{'\b','\t',' '},'MultipleDelimsAsOne',1);
   for i=1:nocol
    block(:,i+1)=blockInCell{i};
   end
  end
  data = [data;block];
 end
 fclose(fid);

end
