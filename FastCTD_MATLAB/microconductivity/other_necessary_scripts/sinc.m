function ans = sinc(x)
% sinc(x) = sin(x)/x.  No pi's.
m=length(x);
ans=zeros(1,m);
for k=1:m,
 if ( abs(x(k)) > .0001),
  ans(k) = sin(x(k))/x(k);
 else
  ans(k) = 1 - x(k)*x(k)/6;
 end
end
