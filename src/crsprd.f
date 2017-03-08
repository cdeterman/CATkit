      subroutine crsprd(n,xb,yb,xc,yc,cp)
c
c     written by Germaine Cornelissen, December 1981
c
c     This subroutine computes cross-products (CP) as-one-goes, using the
c     formula CP(n) = CP(n-1) + (n-1)/n * (xb-xc) * (yb-yc).
c     n indicates the current step being processed
c     xb and yb are the current means of x and y at nth step
c     xc and yc are the current values of x and y : x(n) and y(n).
c
      fn=n
      cp=cp+(fn-1.)/fn*(xb-xc)*(yb-yc)
      xb=((fn-1.)*xb+xc)/fn
      yb=((fn-1.)*yb+yc)/fn
      return
      end
