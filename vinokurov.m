%оценка числовых характеристик
clc;
stolb =9;
N = 500;
xlsdata = xlsread('dann3.xlsx');
data = xlsdata(:,stolb);



data2=data;

minn=-999999;

for i = 1:N
for j=1:N-1
if data2(j+1)<data2(j)
    minn=data2(j+1);
    data2(j+1)=data2(j);
    data2(j)=minn;
end
end
end

   


data2;
if (rem(N,2)==0)
   sr2= data2(ceil(N/2));
   sr1= data2(ceil(N/4));
   sr3= data2(ceil(3*N/4));
   
   mej=sr3-sr1;
   mej=mej*1.5;
   sr11=sr1-mej;
   sr33=sr3+mej;
   
   data3=data;
   
   shet=0;
   for i = 1:N
       
      if (data3(i)<sr11 ||data3(i)>sr33)
          shet=shet+1;
          
          for j=i:N-i
              data3(i)=data3(i+1);
              
          end
          N=N-1;
      end
              
   end
   data2;
   sr33;
   sr11;
   shet;
   N;
   if (shet ~= 0)
       fprintf('В выборке обнаруженны выбросы. В количестве');
       shet
      fprintf(' Они удалены из выборки, её размерность изменена. Теперь её размерность');
    
N
   else
        fprintf('В выборке выбросов не обнаружненно \n');
   end
else
  sr2=(data2(ceil(N/2))+data2(floor(N/2)));
  sr1=data2(ceil(N/4));
  sr3= data2(ceil(3*N/4));
  
   mej=sr3-sr1;
   mej=mej*1.5;
   sr11=sr1-mej;
   sr33=sr3+mej;
   
   data3=data;
   
   shet=0;
   for i = 1:N
       
      if (data3(i)<sr11 ||data3(i)>sr33)
          shet=shet+1;
          
          for j=i:N-i
              data3(i)=data3(i+1);
              
          end
          N=N-1;
      end
              
   end
   data2;
   sr33;
   sr11;
   shet;
   N;
   
   if (shet ~= 0)
       fprintf('В выборке обнаруженны выбросы. В количестве');
       shet
      fprintf(' Они удалены из выборки, её размерность изменена. Теперь её размерность');
       
   else
        fprintf('В выборке выбросов не обнаружненно \n');
   end
   

end

fprintf('____________________________________________________________________________\n');
 


summ=0;
for i = 1:N
 summ =summ+ data3(i);
 end
%summ
sr=0;
sr=summ/N

Dispersiya=0;
summd=0;
for i = 1:N
 summd =summd+data3(i)^2;
end
 Dispersiya=(1/(N-1))*(summd-N*sr*sr);
 Dispersiya
 
 SrKvOtkl=0;
 SrKvOtkl=sqrt(Dispersiya);
 SrKvOtkl
 
 KoefVar=0;
 KoefVar=SrKvOtkl/sr;
 KoefVar
 
 KoefAssim=0;
 m3=0;
  for i = 1:N
 m3=m3+((data3(i)-sr).^3);
 end
 m3=m3*(1/N);
 KoefAssim=m3/(SrKvOtkl.^3);
 KoefAssim
 
 Ekscess=0;
 m4=0;
  for i = 1:N
 m4=m4+((data3(i)-sr).^4);
 end
 m4=m4*(1/N);
 Ekscess=(m4/(SrKvOtkl.^4))-3;
 Ekscess
 
 Median = median(data3);
 Median
 
 Moda=mode(data3);
 Moda
 
 
 
  max=-999999;
 for i = 1:N
     if data(i)>max
 max=data3(i);
     end
 end
max;
 
 
 min=999999;
 for i=1:N
     if data(i)<min
         min=data3(i);
     end
 end
 
 min;
 
 
 
n = fix(sqrt(N)); %на сколько интервалов разбиваем корень из количества

a(1,n) = 0;
mn=0;
karman = (max-min)/n; %длинна интервала
for i = 1:500
    for j = 1:n
       if ((data3(i)>=(min+karman*(j-1)) &&(data3(i)<(min+karman*j))))
         a(j) = a(j)+1;   
       end
    end
    
    if max == data3(i)
        mn = mn+1;
    end
end

if min+karman*j == max
    a(n) = a(n) + mn;
end

 
 
  
 %----------------------------------------------------------
%графический анализ
subplot(2,2,1) 
cdfplot(data3)  
subplot(2,2,2)  
nbins=round(sqrt(N));
h=histfit(data3,nbins);

 
 %----------------------------------------------------------
%предварительная проверка на нормальность
 fprintf('____________________________________________________________________________\n');
%односигмовый интервал
sigma1 = 0;
for i = 1:N
    if (data3(i) > Median - SrKvOtkl) && (data3(i) < Median + SrKvOtkl)
        sigma1 = sigma1 +1;
    end
end


sigma1/N*100;

sigma2 = 0;
for i = 1:N
    if (data3(i) >( Median - 2* SrKvOtkl)) && (data3(i) < (Median + 2* SrKvOtkl))
        sigma2 = sigma2 +1;
    end
end

sigma2/N*100;

sigma3 = 0;
for i = 1:N
    if (data3(i) >( Median - 3* SrKvOtkl)) && (data3(i) < (Median + 3* SrKvOtkl))
        sigma3 = sigma3 +1;
    end
end

sigma3/N*100;

if (sigma1>68 && sigma2>=95 && sigma3 >=99)
    sigma=true;
   fprintf('В односигмовый, двухсигмовый и трёхсигмовый интервал вошли 68 95 и 100 процентов выборки\n')
else
    sigma=false;
fprintf('В односигмовый, двухсигмовый и трёхсигмовый интервал не вошли 68 95 и 100 процентов выборки\n')
end

KoefVar;
if KoefVar < 0.33
    V=true;
   fprintf('Коэффициент вариации V < 0.33\n');
else
    V=false;
   fprintf('Коэффициент вариации V > 0.33\n');
end

KoefAssim;
Ekscess;

if (KoefAssim <=2 ||  KoefAssim>=2 && Ekscess<=2 || Ekscess>=2)
    KoefAssimANDEkscess=true;
   fprintf('Коэффициент ассиметрии и эксцесс близки к нулю \n');
else
    KoefAssimANDEkscess=false;
   fprintf('Коэффициент ассиметрии и эксцесс не близки к нулю \n');
end

sr ;
Median;

if (sr-Median<=1 || sr-Median>=-1)
     rav=true;
   fprintf('Среднее и медиана примерно равны \n');
else
    rav=false;
   fprintf('Среднее и медиана сильно отличаются\n');
end

if(rav==true && KoefAssimANDEkscess==true &&  V==true && sigma==true)
    fprintf('Выборка предварительно нормальна\n');
else
    fprintf('Выборка предварительно не нормальна \n');
end
    
 fprintf('____________________________________________________________________________\n');

Dispersiya;

per=0;
if (N>=20)
    Ckv=0;
    Sx=1/(sqrt(N+1));
    t=0;
     for i = 1:N-1
         per=per+((data3(i+1)-data3(i)).^2);
     end
 
     Ckv=(1/(2*(N-1)))*per;
     t= Ckv/Sx;
     tkr=1.960;
     t;
    if t<tkr  fprintf('Выборка случайна \n');
    else  fprintf('Выборка не случайна \n');
    end
    
     
else
    
     Ckv=0;
    t=0;
     for i = 1:N-1
         per=per+((data3(i+1)-data3(i)).^2);
     end
 
     Ckv=(1/(2*(N-1)))*per;
     t= Ckv/Dispersiya;
     tkr=1.65;
     t;
    if t<tkr  fprintf('Выборка случайна \n');
    else  fprintf('Выборка не случайна \n');
    end
    
end
 fprintf('____________________________________________________________________________\n');



 %тренд
A = 0;
for i=1:N
    for j=i+1:N
        if data3(i) > data3(j)
        A=A+1;   
        end
    end
end
A;

Ua=(N*(N-1))/4;
Ba=(N*(2*N+5)*(N-1)/72);


Za=1.96;

minn=Ua-Za*Ba;
maxx=Ua+Za*Ba;

if (A>minn && A<max)
    fprintf('В выборке тренд не найден \n')
else
    fprintf('В выборке найден тренд \n')
end

fprintf('____________________________________________________________________________\n');

%знакочередавание
A=0;
Aplus=0;
Aminus=0;
for i=1:N-1   
        if ((data3(i+1)<0 &&data3(i)>0) || (data3(i+1)>0 &&data3(i)<0))
            A=A+1;
        end  
end
for i=1:N
    if data3(i)>=0
        Aplus=Aplus+1;
    else
        Aminus=Aminus+1;
    end    
end



Aplus;
Aminus;

A1=(((2*Aplus*Aminus)/N)+1)-Za*(sqrt(2*Aplus*Aminus*((2*Aplus*Aminus)-N)/(N*N*(N-1))));
 
A2=(((2*Aplus*Aminus)/N)+1)+Za*(sqrt((2*Aplus*Aminus*((2*Aplus*Aminus)-N)/N*N*(N-1))));
 


if (A/Za>A1 &&A/Za<A2)
    fprintf('Выборка знакочередующаяяся \n')
else
    if (A/Za*Za>A1 &&A/Za*Za<A2)
        fprintf('Выборка условно знакочередующаяяся \n')
    else
    fprintf('Выборка не знакочередующаяяся \n')
    end
end

    
fprintf('____________________________________________________________________________\n');



Std = std(data);

F = chi2gof(data,'cdf',{@normcdf, sr, Std});  
    if F == 0
        fprintf('\nРаспределение нормальное(хи-квадрат)\n')
    else
        fprintf('\nРаспределение не нормальное(хи-квадрат)\n')
    end

    
    F = chi2gof(data,'cdf',{@unifcdf, min, max});
    if F == 0
        fprintf('Распределение равномерное(хи-квадрат)\n')
    else
        fprintf('Распределение не равномерное(хи-квадрат)\n')
    end

