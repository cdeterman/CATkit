      function phsrd(x,y)
c
c     This function determines the polar angle from rectangular
c     coordinates (x,y).
c
      pi=4.0*atan(1.0)
      dgprd=180./pi
      val=10000.0*abs(x)-abs(y)
      if(val.gt.0.)go to 10
      phsrd=90.0
      go to 20
   10 phsrd=dgprd*atan(y/x)
      if(phsrd.ge.0.)go to 20
      phsrd=phsrd+180.0
   20 if(y)30,40,50
   30 phsrd=phsrd+180.0
      return
   40 phsrd=0.0
      if(x.lt.0.)go to 30
   50 phsrd=phsrd
      return
      end
