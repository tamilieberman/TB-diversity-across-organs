function num=locustag2num(tag)

x=char(tag);
if numel(x) > 1
    num=find(x-0>=48 & x-0<58,1);
else
    num=0;
end

end