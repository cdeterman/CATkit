      function gamlog(z)
c
c     This function computes the log of the gamma function of z. Computations 
c     are based on a polynomial approximation given by M. Abramowitz and
c     A. Stegun (Handbook of Mathematical Functions, Dover, 1965;
c     formula 6.1.36).
c

      real  :: z
      
      if(z.le.0.)stop
      iz=ifix(z)
      x=z-float(iz)
      gamlog=1.+x*(-.577191652+x*(.988205891+x*(-.897056937+x*(.918206857
     *+x*(-.756704078+x*(.482199394+x*(-.193527818+x*.035868343)))))))
      if(z.gt.1.)go to 10
      gamlog=gamlog/z
      gamlog=alog(gamlog)
      return
   10 i=2
      gamlog=alog(gamlog)
   20 if(i.gt.iz)return
      x=x+1.
      gamlog=alog(x)+gamlog
      i=i+1
      go to 20
      end
