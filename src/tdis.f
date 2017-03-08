      function tdis(z,ndf)
c
c     This function calls 'beta' to compute the beta function.
c
      if(ndf.le.0)stop
      df=ndf
      zz=1./(1.+z*z/df)
      df=.5*df
      tdis=1.-.5*beta(zz,df,.5)
      if(z.lt.0.)tdis=1.0-tdis
      return
      end
