      subroutine phase(arg,c33,c23,acr,angle)
c
c     written by Germaine Cornelissen, December 1981
c
c     This subroutine determines the limits of the (1-alpha)%confidence 
c     interval for the acrophase.
c
      pi=4.0*atan(1.0)
      phi=atan(arg)
      xt=c33*cos(phi)+c23*sin(phi)
      if(xt.lt.0.)phi=phi+pi
      phi=phi*180./pi
      angle=phi+acr
      if(angle.gt.0.)angle=angle-360.
      if(angle.lt.-360.)angle=angle+360.
      return
      end
