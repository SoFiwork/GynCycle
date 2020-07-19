function S =poissonproc(lambda,tspan)

tb=tspan(1);    %start time
te=tspan(2);  %end time

%alternative method:
N=poissrnd(lambda*(te-tb));
U=rand([N,1]);
Usort=sort(U);
S=tb+(te-tb)*Usort;
%length(S)

end
