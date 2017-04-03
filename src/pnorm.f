      function pnorm(z)
c
c     This function computes the probability associated with the standard
c     normal distribution. The polynomial approximation used is that given
c     by M. Abramowitz and A. Stegun (Handbook of Mathematical Functions,
c     Dover, 1965; formula 26.2.19).
c
      real  :: z
      
      x=abs(z)
      x=(((((.0000053830*x+.0000488906)*x+.0000380036)*x
     *+.0032776263)*x+.0211410061)*x+.0498673470)*x+1.0
      x=x*x
      x=x*x
      x=x*x
      x=x*x
      x=1.0-0.5/x
      if(z.lt.0.)x=1.0-x
      pnorm=x
      return
      end