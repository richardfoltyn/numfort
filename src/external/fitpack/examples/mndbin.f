cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc              mndbin : dblint test program                          cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 fac,facx,aint,exint,dblint,xb,xe,yb,ye
      integer i,j,kx,kx1,ky,ky1,mx,my,m0,m1,m2,m3,nc,nkx1,nky1,nx,ny
      real*8 tx(15),ty(15),c(100),x(6),y(6),wrk(50)
c  we set up the end points of the integration domains.
      mx = 6
      my = 6
      do 10 i=1,6
      x(i) = (i-1)*0.2
      y(i) = x(i)
  10  continue
c  loop for different spline degrees with respect to the x-variable
      do 300 kx=1,5,2
c  the knots in the x-direction
        tx(kx+2) = 0.4
        tx(kx+3) = 0.7
        tx(kx+4) = 0.9
        kx1 = kx+1
        nx = 3+2*kx1
        j = nx
        do 20 i=1,kx1
          tx(i) = 0.
          tx(j) = 1.
          j = j-1
  20    continue
c  loop for different spline degrees with respect to the y-variable
      do 200 ky=2,3
c  the knots in the y-direction
        ty(ky+2) = 0.3
        ty(ky+3) = 0.8
        ky1 = ky+1
        ny = 2+2*ky1
        j = ny
        do 30 i=1,ky1
          ty(i) = 0.
          ty(j) = 1.
          j = j-1
  30    continue
c  we generate the b-spline coefficients for the test function x*y
        nkx1 = nx-kx1
        nky1 = ny-ky1
        do 40 i=1,nky1
          c(i) = 0.
  40    continue
        do 50 i=2,nkx1
          c((i-1)*nky1+1) = 0.
  50    continue
        fac = kx*ky
        m0 = 1
        do 70 i=2,nkx1
          m1 = m0+nky1
          facx = (tx(i+kx)-tx(i))/fac
          do 60 j=2,nky1
            m2 = m0+1
            m3 = m1+1
            c(m3) = c(m1)+c(m2)-c(m0)+facx*(ty(j+ky)-ty(j))
            m0 = m0+1
            m1 = m1+1
  60      continue
          m0 = m0+1
  70    continue
c  printing of the spline information
        write(6,900) kx,ky
        write(6,910)
        write(6,920) (tx(i),i=1,nx)
        write(6,930)
        write(6,920) (ty(i),i=1,ny)
        nc = nkx1*nky1
        write(6,940)
        write(6,950) (c(i),i=1,nc)
c  integration of the tensor product spline
        write(6,960)
        j = 6
        do 100 i=1,3
          xb = x(i)
          yb = xb
          xe = x(j)
          ye = xe
          j = j-1
          aint = dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk)
          exint = (xe-xb)*(xe+xb)*(ye-yb)*(ye+yb)*0.25
          write(6,970) xb,xe,yb,ye,aint,exint
 100    continue
 200    continue
 300  continue
      stop
c  format statements.
 900  format(33h0tensor product spline of degrees,2i3)
 910  format(1x,40hposition of the knots in the x-direction)
 920  format(1x,15f5.1)
 930  format(1x,40hposition of the knots in the y-direction)
 940  format(23h b-spline coefficients )
 950  format(1x,8f9.4)
 960  format(1h0,5x,2hxb,5x,2hxe,5x,2hyb,5x,2hye,6x,6hdblint,7x,5hexint)
 970  format(1x,4f7.1,2f12.5)
      end
