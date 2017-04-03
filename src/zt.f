      function ztdis(p,ndf)
c
c     This function computes the value of the t-distribution given the
c     probability level and the number of degrees of freedom. Other
c     functions required are 'tdis', 'beta' and 'gamlog'.
c

      real, intent(in)       ::  ndf
      real                   ::  p
      
      n=0
      if(ndf.le.0)stop
      if(p.lt.0..or.p.gt.1.)stop
      zz=0.
      del=p
      if(p.ge.0.5)go to 10
      p=1.-p
      del=p
      zz=1.
   10 if(ndf.ge.4)go to 20
      z3=ndf
      del=(12.*p*p)/(z3*z3)
   20 z3=0.
      p3=0.5
   30 z1=z3
      p1=p3
      z3=z1+del
      p3=tdis(z3,ndf)
      if(p3-p)30,90,40
   40 z2=z3
      p2=p3
   50 dz=(z2-z1)*(p-p1)/(p2-p1)
      z3=z1+dz
      p3=tdis(z3,ndf)
      if(p3-p)60,90,70
   60 d=p-p3
      z1=z3
      p1=p3
      go to 80
   70 d=p3-p
      z2=z3
      p2=p3
   80 if(d.lt.0.00005)go to 90
      n=n+1
      if(n.le.60)go to 50
   90 if(zz.lt.1.)go to 100
      z3=-z3
      p=1.-p
  100 ztdis=z3
      return
      end
