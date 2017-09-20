function [ Data_clean ] = cleanData( Data)
%this function eliminates noise zero countour from data
Ndata=length(Data);
Data_clean=zeros(Ndata,1)
for i=1:Ndata
    if(abs(Data(i,1))<10^-2)
        Data_clean(i,1)=0.;
    else
        Data_clean(i,1)=Data(i,1);
    end
end

