!
! A routine to erase the corresponding output when add=0
!
! Published in [ Jon F. Claerbout, 1992, Earth Soundings Analysis:
!                Processing versus Inversion, Blackwell Scientific Publication ]
!
! Author : Jon F. Claerbout, claerbout@stanford.edu, 1992
subroutine adjnull( adj, add, x, nx,  y, ny )
integer ix, iy,     adj, add,    nx,     ny
real                          x( nx), y( ny )
if ( add .eq. 0 ) then
  if ( adj .eq. 0 ) then
    do iy= 1, ny
      y(iy) = 0.
    end do
  else
    do ix= 1, nx
      x(ix) = 0. 
    end do
  end if
end if
return
end
!
! A function to generate a random number
!
! Published in [ Jon F. Claerbout, 1992, Earth Soundings Analysis:
!                Processing versus Inversion, Blackwell Scientific Publication ]
!
! Author : Jon F. Claerbout, claerbout@stanford.edu, 1992
real function rand01( iseed)
integer ia, im,       iseed
parameter(ia = 727,im = 524287)
iseed = mod(iseed*ia,im)
rand01 = (float(iseed) - 0.5)/float(im - 1)
return
end
!
! A step of conjugate-gradient descent algorithm.
!
! Published in [ Jon F. Claerbout, 1992, Earth Soundings Analysis:
!                Processing versus Inversion, Blackwell Scientific Publication ]
!
! Author : Jon F. Claerbout, claerbout@stanford.edu, 1992
subroutine cgstep( iter,   n, x, g, s,   m, rr, gg, ss)
integer i, n, m, iter
real      x(n), rr(m),  g(n), gg(m),  s(n), ss(m)
! solution, residual
! gradient, conjugate gradient
! step,     conjugate step
real dot, sds, gdg, gds, determ, gdr, sdr, alfa, beta
if ( iter .eq. 0 ) then
  do i= 1, n
    s(i) = 0.
  end do 
  do i= 1, m
    ss(i) = 0.
  end do 
  if ( dot(m,gg,gg).eq.0 ) then
    call erexit('cgstep: grad vanishes identically')
  end if
  alfa =  - dot(m,gg,rr) / dot(m,gg,gg)
  beta = 0.
else
! search plane by solving 2-by-2
  gdg = dot(m,gg,gg)      !  G . (R - G*alfa - S*beta) = 0
  sds = dot(m,ss,ss)      !  S . (R - G*alfa - S*beta) = 0
  gds = dot(m,gg,ss)
  determ = gdg * sds - gds * gds + 1.e-15
  gdr = - dot(m,gg,rr)
  sdr = - dot(m,ss,rr)
  alfa = ( sds * gdr - gds * sdr ) / determ
  beta = (-gds * gdr + gdg * sdr ) / determ
end if 
do i= 1, n                                      ! s = model step
  s(i)  = alfa * g(i)  + beta * s(i)
end do 
do i= 1, m                                      ! ss = conjugate
  ss(i) = alfa * gg(i) + beta * ss(i)
end do 
do i= 1, n                                      ! update solution
  x(i)  =  x(i) +  s(i)
end do 
do i= 1, m                                      ! update residual
  rr(i) = rr(i) + ss(i)
end do 
return
end
real function dot( n, x, y )
integer i, n
real val, x(n), y(n)
val = 0.
do i= 1, n
  val = val + x(i) * y(i)
end do 
dot = val
return
end
!
!
subroutine copy( n, xx, yy)
integer i, n
real xx(n), yy(n)
do i= 1, n
  yy(i) = xx(i)
end do 
return
end
!
subroutine zero(n, xx)
integer i, n
real xx(n)
do i= 1, n
  xx(i) = 0.
end do 
return
end  
