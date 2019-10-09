! -------------------------------------------------------------------
!
! Copyright (c) 2009 by the Society of Exploration Geophysicists.
! For more information, go to http://software.seg.org/2009/0003 .
! You must read and accept usage terms at:
! http://software.seg.org/disclaimer.txt before use.
!
! -------------------------------------------------------------------
!
!       2D RTM algorithm Dot-product test
!
!       test_dot2d.x nodoc
!%
program zyxabc
implicit none
real, allocatable, dimension (:,:) :: data
real, allocatable, dimension (:,:) :: modl
real, allocatable, dimension (:,:) :: dcop
real, allocatable, dimension (:,:) :: mcop
real, allocatable, dimension (:,:) :: dhat
real, allocatable, dimension (:,:) :: mhat
real, allocatable, dimension (:,:) :: v
integer n1,n2,n3
integer getch,putch
call initpar()
call doc('./../Src/test_dot2d.f90')
  if (0>= getch('n1','d',n1 )) then
    n1 = 120
  end if
  if (0.ne.putch('Read  from param: #n1 ','d',n1)) then
    call erexit('Trouble writing n1 to history')
  end if
  if (0>= getch('n2','d',n2 )) then
    n2 = 180
  end if
  if (0.ne.putch('Read  from param: #n2 ','d',n2)) then
    call erexit('Trouble writing n2 to history')
  end if
  if (0>= getch('n3','d',n3 )) then
    n3 = 180
  end if
  if (0.ne.putch('Read  from param: #n3 ','d',n3)) then
    call erexit('Trouble writing n3 to history')
  end if
allocate (data(n1,n2),modl(n3,n2))
allocate (dcop(n1,n2),mcop(n3,n2))
allocate (dhat(n1,n2),mhat(n3,n2))
allocate (v(n3,n2))
call doit ( n1,n2,n3,data,modl,dcop,mcop,dhat,mhat,v)
end program zyxabc 
!
subroutine doit( n1,n2,n3,data,modl,dcop,mcop,dhat,mhat,v)
integer n1,n2,n3,i1,i2,i3, iseed
real rand01, dt,dx
real data(n1,n2),modl(n3,n2)
real dcop(n1,n2),mcop(n3,n2)
real dhat(n1,n2),mhat(n3,n2)
real v(n3,n2)
real dotd, dotm
iseed = 2003
!
! fill d=data(i1,i2) with randum numbers and save a copy to dcop(i1,i2)
!
do i1=1,n1
  do i2=1,n2  
    data(i1,i2) = rand01(iseed)
    dcop(i1,i2) = data(i1,i2)
  end do 
!
! fill m=modl(i3,i2) with randum numbers and save a copy to mcop(i3,i2)
!
end do 
do i2=1,n2
  do i3=1,n3  
    modl(i3,i2) = rand01(iseed)
    mcop(i3,i2) = modl(i3,i2)
  end do 
!
! prepare randum velocity field
!
end do 
do i2=1,n2
  do i3=1,n3
    v(i3,i2) = 2000.*rand01(iseed) + 1000.
  end do
end do 
dt=.004
dx=25.
!
!  d' = L m
!
call rt2d( 0, dt,dx, n1,n2,n3, modl,dhat,v)
call copy(n3*n2,mcop,modl)
!
!  m' = L'd
!
call rt2d( 1, dt,dx, n1,n2,n3, mhat,data,v)
call copy(n1*n2,dcop,data)
dotd = 0.
dotm = 0.
!
!  d' (dot) d
!
do i1=1,n1
  do i2=1,n2
    dotd = dotd + dhat(i1,i2) * data(i1,i2)
  end do 
!
!  m (dot) m'
!
end do 
do i2=1,n2
  do i3=1,n3
    dotm = dotm + modl(i3,i2)*mhat(i3,i2)
  end do
end do 
write(0, 10) "d'd =", dotd
write(0, 10) "m m'=", dotm
10 format( a, f15.8)
return
end  
