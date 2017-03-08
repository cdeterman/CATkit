      function beta(z,a,b)
c
c     This function computes the beta function. The function calls 'gamlog'
c     and 'pnorm'.
c
      ike=1
      if(z.lt.0.)stop
      if(z.gt.0.)go to 10
      beta=0.
      return
   10 if(z.gt.1.)stop
      if(z.lt.1.)go to 20
      beta=1.
      return
   20 if(a.le.0..or.b.le.0.)stop
      if(z.gt.0.84)go to 30
      zz=z
      aa=a
      bb=b
      go to 40
   30 zz=1.-z
      aa=b
      bb=a
      ike=-1
   40 sum=1.
      p=1.
      w1=aa+bb
      w2=aa+1.
      do 50 i=1,6
         p=p*zz*w1/w2
         sum=sum+p
         w1=w1+1.
         w2=w2+1.
   50 continue
      sum=alog(sum)+gamlog(aa+bb)-(gamlog(aa+1.)+gamlog(bb))
      sum=exp(sum)
      sum=sum*(zz**aa)*((1.-zz)**bb)
      aa=w2
      w1=(bb*zz)**(.33333333333)
      w2=(aa*(1.-zz))**(.33333333333)
      aa=1./aa
      bb=1./bb
      zz=sqrt(w1*w1*bb+w2*w2*aa)
      zz=3.*(w1*(1.-.11111111111*bb)-w2*(1.-.11111111111*aa))/zz
      sum=sum+pnorm(zz)
      if(ike.lt.0)sum=1.-sum
      beta=sum
      return
      end
