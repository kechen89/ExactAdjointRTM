! -------------------------------------------------------------------
!
! Copyright (c) 2009 by the Society of Exploration Geophysicists.
! For more information, go to http://software.seg.org/2009/0003 .
! You must read and accept usage terms at:
! http://software.seg.org/disclaimer.txt before use.
!
! -------------------------------------------------------------------
!
!	1D RTM algorithm Dot-product test
!
!	usage: test_dot1d.x nodoc
!%
program zyxabc
implicit none
real, allocatable, dimension (:) :: data
real, allocatable, dimension (:) :: modl
real, allocatable, dimension (:) :: dcop
real, allocatable, dimension (:) :: mcop
real, allocatable, dimension (:) :: dhat
real, allocatable, dimension (:) :: mhat
real, allocatable, dimension (:) :: vel
integer n1,n3
integer getch,putch
call initpar()
call doc('./../Src/test_dot1d.f90')
  if (0>= getch('n1','d',n1 )) then
    n1 = 100
  end if
  if (0.ne.putch('Read  from param: #n1 ','d',n1)) then
    call erexit('Trouble writing n1 to history')
  end if
  if (0>= getch('n3','d',n3 )) then
    n3 = 100
  end if
  if (0.ne.putch('Read  from param: #n3 ','d',n3)) then
    call erexit('Trouble writing n3 to history')
  end if
allocate (data(n1),modl(n3))
allocate (dcop(n1),mcop(n3))
allocate (dhat(n1),mhat(n3))
allocate (vel(n3))
call doit ( n1,n3,data,modl,dcop,mcop,dhat,mhat, vel)
end program zyxabc 
!
subroutine doit( n1,n3,data,modl,dcop,mcop,dhat,mhat, vel)
integer n1,n3,i1,i3, iseed
real rand01
real data(n1),modl(n3)
real dcop(n1),mcop(n3)
real dhat(n1),mhat(n3)
real vel(n3)
real dotd, dotm, dt,dz
iseed = 1992
!
! fill d = data(i1) with randum numbers and save a copy to dcop(i1)
!
do i1=1,n1  
  data(i1) = rand01(iseed)
  dcop(i1) = data(i1)
end do 
!
! fill m=modl(i3) with randum numbers and save a copy to mcop(i3)
!
do i3=1,n3  
  modl(i3) = rand01(iseed)
  mcop(i3) = modl(i3)
end do 
!
! prepare randum velocity field
!
do i3=1,n3  
  vel(i3) = rand01(iseed)*3000.+1000.
end do 
dt=.004
dz=25.
!
!  d' = L m
!
call rt1d( 0, dt,dz, n1,n3, modl,dhat, vel,0 )
call copy(n3,mcop,modl)
!
!  m' = L'd
!
call rt1d( 1, dt,dz, n1,n3, mhat,data, vel,0 )
call copy(n1,dcop,data)
dotd = 0.
dotm = 0.
!
! d' (dot) d
!
do i1=1,n1
  dotd = dotd + dhat(i1) * data(i1)
end do 
!
! m (dot) m'
!
do i3=1,n3
  dotm = dotm + modl(i3)*mhat(i3)
end do 
write(0, 10) "d'd =",dotd
write(0, 10) "m m'=",dotm
10 format(a, f13.8)
return
end  
