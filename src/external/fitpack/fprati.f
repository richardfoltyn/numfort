      pure function fprati(p1,f1,p2,f2,p3,f3)
c  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
c  gives the value of p such that the rational interpolating function
c  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
c  ..
c  ..scalar arguments..
c      real*8 p1,f1,p2,f2,p3,f3
      real*8, intent(in)  :: p1
      real*8, intent(in)  :: f1
      real*8, intent(in)  :: p2
      real*8, intent(in)  :: f2
      real*8, intent(in)  :: p3
      real*8, intent(in)  :: f3

      real*8 fprati
c  ..local scalars..
      real*8 h1,h2,h3,p
c  ..
      if(p3.gt.0.) go to 10
c  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 40
c  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
      ! RF: No idea why this is here, most call sites do not use the
      ! variables p1,f1,etc. after call to fprati, so reordering
      ! these would have had no effect.
c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
c  20  if(f2.lt.0.) go to 30
c      p1 = p2
c      f1 = f2
c      go to 40
c  30  p3 = p2
c      f3 = f2
  40  fprati = p
      return
      end
