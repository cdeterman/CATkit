      subroutine module(a,b,r)
c
c     written by Germaine Cornelissen, December 1981
c
c     This subroutine computes the module of a vector with rectangular
c     projections (a,b).
c
      c=amax1(abs(a),abs(b))
      d=amin1(abs(a),abs(b))
      q=d/c
      r=c*sqrt(1.+q*q)
      return
      end
