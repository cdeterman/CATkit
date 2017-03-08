      subroutine pnorm(z)
      implicit none
      
      double precision, intent(inout)   ::  z
c
c     This function computes the probability associated with the standard
c     normal distribution. The polynomial approzimation used is that given
c     by M. Abramowitz and A. Stegun (Handbook of Mathematical Functions,
c     Dover, 1965; formula 26.2.19).
c
      z=abs(z)
      z=(((((.0000053830*z+.0000488906)*z+.0000380036)*z
     *+.0032776263)*z+.0211410061)*z+.0498673470)*z+1.0
      z=z*z
      z=z*z
      z=z*z
      z=z*z
      z=1.0-0.5/z
      if(z.lt.0.)z=1.0-z
c     pnorm=z
      return
      end
