c     written by Germaine Cornelissen, December 1980
c     last revised January 1987
c     revised by Charles Determan for R interface, March 2017
c
c     This program computes point and (1-alpha)%interval estimates of 
c     rhythm population parameters and tests for statistical significance.
c     The P-value given corresponds to a zero-amplitude test. 
c     Confidence intervals for amplitude and acrophase alone are based
c     on approximations valid only when P is well below alpha.
c
c     The program calls subroutines crsprd, module and phase:
c    -crsprd computes cross-products as-one-goes;
c    -module computes amplitude from rectangular projections;
c    -phase determines the limits of the (1-alpha)%confidence interval  
c     for the acrophase.
c
c     The program also calls library routines phsrd, setfil and p2p:
c    -phsrd computes the acrophase from beta and gama;
c    -setfil and p2p are machine-dependent plotting routines and should be
c     replaced by user's own software. All instructions relative to machine-
c     dependent plotting software are braketed by c$$$ cards.
c
c     Labels 1000's are usually reserved to read and 2000's to write.
c     Input files are 5 (standard) and iu (e.g., iu=10) if setup and data
c     files are separate. Output files are 6 (standard) and 7 (if copy of
c     rhythm parameters from input is requested). Files 14, 15 and 19 are
c     used in relation to the cosinor plot.
c
c     sig(i,j) are the elements of the variance-covariance matrix; c22,c23
c     and c33 are elements of that matrix after rotation, used to compute
c     confidence intervals for amplitude and acrophase according to Fieller's
c     method; rytpar is an array containing information for the cosinor plot.
c

      SUBROUTINE bmpopsub(rmatrix, M, N, sig, rytpar,
     $alpha, k, perd)
      
      IMPLICIT NONE
      
      integer, intent(in), value      ::  k, M, N
      real, intent(in)      ::  rmatrix(M,N)
      real, intent(in), value   :: perd
      real, intent(inout)   ::  alpha
      real, intent(inout)   ::  sig(3,3)
      real, intent(inout)   ::  rytpar(6,13)
c     real, intent(inout)   ::  ames(4), 
c    &      btsum(4), gmsum(4)
      real  :: ames(4), btsum(4), gmsum(4)
     
      real :: ia1, ia2, cim, ciamp, c, c22, c33, c23,
     &      den, fdd, d, btmes, btgam, bt2sum, beta, b, arg, 
     &      an1, an2, phsrd, ames2, acr, a, a1, a2, p, pi, phi1, 
     &      phi2, ndf, ia, ipr, amp, x1, t, sg, r, prs, pr, pp,
     &      plev, ztdis, ymes, x2, gm2sum, games, gama,
     &      fval, fk, f, df,dd
      integer :: i
      
      pi=4.0*atan(1.0)
      if(alpha <= 0) alpha=0.05
      plev=100.*(1.-alpha)
c
c     do 132 kk=1,npop
c
      prs=0.
      do 50 i=1,4
         ames(i)=0.
         btsum(i)=0.
         gmsum(i)=0.
   50 continue
      ames2=0.
      bt2sum=0.
      gm2sum=0.
      btmes=0.
      games=0.
      btgam=0.
c
c*** computation of population parameter estimates and 
c***             of variance-covariance matrix
c
c     if(iprt.eq.1)write(7,1030)titl2
      do 60 i=1,k
c        read(iu,frmt)pr,ymes,amp,acr
c        if(iprt.eq.1)write(7,frmt)pr,ymes,amp,acr
         pr = rmatrix(i,1)
         ymes = rmatrix(i,2)
         amp = rmatrix(i,3)
         acr = rmatrix(i,4)
         
         acr=acr*pi/180.
         beta= amp*cos(acr)
         gama=-amp*sin(acr)
         prs=prs+pr
         call crsprd(i,ames(1) ,ames(2) ,ymes,ymes,ames2 )
         call crsprd(i,btsum(1),btsum(2),beta,beta,bt2sum)
         call crsprd(i,gmsum(1),gmsum(2),gama,gama,gm2sum)
         call crsprd(i,ames(3) ,btsum(4),ymes,beta,btmes )
         call crsprd(i,btsum(3),gmsum(4),beta,gama,btgam )
         call crsprd(i,gmsum(3),ames(4) ,gama,ymes,games )
   60 continue
      fk=k
      pr=prs/fk
	ipr=pr + 0.5
      ymes=ames(1)
      beta=btsum(1)
      gama=gmsum(1)
      call module(beta,gama,amp)
      acr=-phsrd(beta,gama)
      df=fk-1.
      sig(1,1)= ames2/df
      sig(1,2)= btmes/df
      sig(1,3)= games/df
      sig(2,2)=bt2sum/df
      sig(2,3)= btgam/df
      sig(3,3)=gm2sum/df
      sig(2,1)=sig(1,2)
      sig(3,1)=sig(1,3)
      sig(3,2)=sig(2,3)
c
c*** computation of confidence intervals for the population parameters
c
      ndf=k-1
      sg=1.-alpha/2.
      t=ztdis(sg,ndf)
      cim=t*sqrt(sig(1,1)/fk)
      c22=sig(2,2)*beta*beta+2.*sig(2,3)*beta*gama+sig(3,3)*gama*gama
      c23=(sig(3,3)-sig(2,2))*beta*gama+sig(2,3)*(beta*beta-gama*gama)
      c33=sig(2,2)*gama*gama-2.*sig(2,3)*beta*gama+sig(3,3)*beta*beta
      den=fk*amp*amp
      c22=c22/den
      c23=c23/den
      c33=c33/den
      if(c22.gt.0.)go to 70
      ciamp=-0.1
      go to 80
   70 ciamp=t*sqrt(c22)
   80 den=amp*amp-c22*t*t
      an2=amp*amp-(c22*c33-c23*c23)*t*t/c33
      if(an2.lt.0.)go to 90
      an2=sqrt(an2)
      an2=t*sqrt(c33)*an2
      an1=c23*t*t
      arg=(an1+an2)/den
      call phase(arg,c33,c23,acr,phi1)
      arg=(an1-an2)/den
      call phase(arg,c33,c23,acr,phi2)
      go to 100
   90 phi1=0.
      phi2=0.
  100 if(ciamp.gt.0.)go to 110
      a1=0.
      a2=0.
      go to 120
  110 a1=amp-ciamp
      a2=amp+ciamp
  120 ia=acr-0.5
      ia1=phi1-0.5
      ia2=phi2-0.5
c
c*** computation of the (1-alpha)%confidence ellipse
c
      f=((1./alpha)**(2./(fk-2.))-1.)*(fk-2.)/2.
      a=1./sig(2,2)
      b=-2.*sig(2,3)/(sig(2,2)*sig(3,3))
      c=1./sig(3,3)
      r=sig(2,3)/sqrt(sig(2,2)*sig(3,3))
      fval=fk*(fk-2.)/(2.*(fk-1.)*(1.-r*r))*(beta*beta/sig(2,2)
     *-2.*r*beta*gama/sqrt(sig(2,2)*sig(3,3))+gama*gama/sig(3,3))
      p=(1.+2.*fval/(fk-2.))**(-(fk-2.)/2.)
      pp=x1
      if(p.gt.0.0005)go to 124
      p=0.001
      pp=x2
  124 if(p.le.alpha)go to 125
      a1=0.
      a2=0.
      ia1=0
      ia2=0
c 125 write(6,2040)kk,titl2,k,ipr,pp,p,ymes,cim,amp,a1,a2,ia,ia1,ia2,
c     *perd
  125 d=2.*(fk-1.)*(1.-r*r)/(fk*(fk-2.))*f
      dd=(fk-1.)*(1.-r*r)/(fk*(fk-2.))*t*t
c
c     Equation of (1-alpha)%confidence ellipse is 
c     a*(x-beta)**2 + b*(x-beta)*(y-gama) + c*(y-gama)**2 = d (or dd)
c     Note that two ellipse equations are given. The first one (using d)
c     provides the correct (1-alpha)%joint confidence region for (A,phi), 
c     but conservative intervals for separate amplitude and acrophase when 
c     taking distances and tangents from the pole to the ellipse. By using
c     dd instead, approximate (1-alpha)%confidence intervals for amplitude
c     and acrophase considered alone are obtained.
c
      do 130 i=1,2
         rytpar(i, 1)=perd
         rytpar(i, 2)=pr
         rytpar(i, 3)=p
         rytpar(i, 4)=amp
         rytpar(i, 5)=ia
         rytpar(i, 6)=a1
         rytpar(i, 7)=a2
         rytpar(i,11)=ymes
         rytpar(i,12)=cim
         rytpar(i,13)=fk
  130 continue
c     rytpar(1, 8)=a/d
      rytpar(1, 8)=ia
c     rytpar(1, 9)=b/d
      rytpar(1, 9)=ia1
c     rytpar(1,10)=c/d
      rytpar(1,10)=ia2
      rytpar(2, 8)=a/dd
      rytpar(2, 9)=b/dd
      rytpar(2,10)=c/dd
c
c 132 continue
c
      end
