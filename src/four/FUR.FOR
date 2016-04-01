      subroutine fftsh (m,x1,x,y)
c     implicit double precision(a-h,o-z)
      dimension x(0:1023),amp(0:1023),ph(0:1023),y(0:1023),x1(1024)
      integer o,e,f
      n=2**m
      do 1 i=1,n
    1 x(i-1)=x1(i)            
          pi=acos(-1.)
      do 59 i=0,n-1
  59  y(i)=0.
           do 17 l=1,m
           e=2**(m+1-l)
           f=e/2
           af=f
           u=1
           v=0
           z=pi/af
           c=cos(z)
           s=-sin(z)
               do 19 j=1,f
                   do 23 i=j,n,e
                   o=i+f-1
                   p=x(i-1)+x(o)
                   q=y(i-1)+y(o)
                   r=x(i-1)-x(o)
                   t=y(i-1)-y(o)
                   x(o)=r*u-t*v
                   y(o)=t*u+r*v
                   x(i-1)=p
                   y(i-1)=q
 23                continue
              w=u*c-v*s
              v=v*c+u*s
              u=w
 19           continue
 17       continue
                     j=1
                     do 47 i=1,n-1
                     if(i.ge.j) goto 235
                     j1=j-1
                     i1=i-1
                     p=x(j1)
                     q=y(j1)
                     x(j1)=x(i1)
                     y(j1)=y(i1)
                     x(i1)=p
                     y(i1)=q
 235                 k=n/2
 240                 continue
                           if(k.ge.j) goto 260
                           j=j-k
                           k=k/2
                           goto 240
 260                j=j+k
 47                 continue
            do 325 k=0,n-1
            a=sqrt(x(k)**2+y(k)**2)
            q=0
                 amp(k)=2*a/n
                 if(abs(a).lt.1.e-7)   goto 324
                 q=acos(x(k)/a)
                    if(y(k).lt.0)  q=-q
 324             ph(k)=q
c       write(6,*) x(k),y(k)
 325        continue
            amp(0)=amp(0)/2
           return
                  end
