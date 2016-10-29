function diversity_from_pileup(fname_in,fname_out,RemoveEnds,ChrStarts,GenomeLength)

tic ;


data=zeros(9,GenomeLength); %[A T C G a t c g D]
% D is number of reads supporting indels

fid=fopen(fname_in);

nts='ATCGatcg';
indels='+-';

line = fgetl(fid);
while ischar(line)
    temp=zeros(9,1); k=0;
    lt=textscan(line,'%s','\t');
    l=lt{1};
    chr=l{1};
    position= ChrStarts((str2double(chr(28))+1))+str2double(l{2});
    ref=strfind(nts,l{3});
    %coverage=str2double(l{4});     quality=l{6};
    calls=l{5};
    while k < length(calls)
        k=k+1;
        c=calls(k);
        nt=strfind(nts,c);
        id=strfind(indels,c);
        if c=='^' %end of read, next character is a mapping quality character
            k=k+1+RemoveEnds;
        elseif c=='$'
            k=k+RemoveEnds; %start of read, skip next character if RemoveEnds is on
        elseif nt
            temp(nt)=temp(nt)+1;
        elseif c=='.'
            temp(ref)=temp(ref)+1;
        elseif c==','
            temp(ref+4)=temp(ref+4)+1;
        elseif id
            k=k+2+str2double(calls(k+1));
            temp(9)=temp(9)+1;
        end
    end
    data(:,position)=temp;
    disp(line)
    disp(temp)
    line = fgetl(fid);
    
end

fclose(fid);

save(fname_out,'data');
            
            
            
