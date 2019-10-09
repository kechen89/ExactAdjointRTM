! -------------------------------------------------------------------
!
! Copyright (c) 2009 by the Society of Exploration Geophysicists.
! For more information, go to http://software.seg.org/2009/0003 .
! You must read and accept usage terms at:
! http://software.seg.org/disclaimer.txt before use.
!
! -------------------------------------------------------------------
!
!  2-D Time-extrapolation modeling and Reverse-time Migration (RTM)
!  using subroutine rt2d_conv() where RTM with conventional implementation
!
!  Input
!       adj=[0/1] : operation type indicator ( 0 for forward / 1 for adjoint )
!       in.H : SEPlib header file for 2D wavefield
!              should contains n1:[nz/nt], d1:[dz/dt], n2:nx, d2:dx
!       vel=vel.H : SEPlib header file for the 2D velocity function
!                   should contains n1:nz, d1:dz, n2:nx, d2:dx
!
!  Output
!       out.H : data [when adj=0] or image [when adj=1]
!  Usage
!
!       test_rt2d_conv.x < in.H adj=[0/1] vel=vel.H > out.H
!
!%----------------
program zyxabc
implicit none
real, allocatable, dimension (:,:) :: modl
real, allocatable, dimension (:,:) :: v
real, allocatable, dimension (:,:) :: data
integer n1,n2,n3, nt,nx,nz, adj
integer getch,hetch,putch,auxpar
call initpar()
call doc('./../Src/test_rt2d_conv.f90')
  if (0>= getch('adj','d',adj )) then
    adj = 0
  end if
  if (0.ne.putch('Read  from param: #adj ','d',adj)) then
    call erexit('Trouble writing adj to history')
  end if
if (adj.eq.0) then
    if (0>= hetch('n1','d',nz )) then
      call erexit('Trouble reading  n1 from history  ')
    end if
    if (0>= hetch('n2','d',nx )) then
      call erexit('Trouble reading  n2 from history  ')
    end if
    if (0>= getch('nt','d',nt )) then
      call erexit('Trouble reading  nt from param  ')
    end if
    if (0.ne.putch('Read  from param: #nt ','d',nt)) then
      call erexit('Trouble writing nt to history')
    end if
else
    if (0>= hetch('n1','d',nt )) then
      call erexit('Trouble reading  n1 from history  ')
    end if
    if (0>= hetch('n2','d',nx )) then
      call erexit('Trouble reading  n2 from history  ')
    end if
    if (0>= auxpar('n1','d',nz ,'vel')) then
      call erexit('Trouble reading  n1 from aux vel ')
    end if
    if (0.ne.putch('Read  from aux: vel_n1 ','d',nz)) then
      call erexit('Trouble writing nz to history')
    end if
end if 
allocate (modl(nz,nx))
allocate (v(nz,nx))
allocate (data(nt,nx))
call submain ( nt,nx,nz, modl,data, v, adj)
end program zyxabc 
!
subroutine submain( nt,nx,nz, modl,data, v, adj)
integer n1,n2,n3, nt,nx,nz, i,j, adj
real    d1,d2,d3, dt,dx,dz, o1,o2,o3 
real    modl(nz,nx)
real    v(nz,nx)
real    data(nt,nx)
integer getch,hetch,putch,auxpar
if (adj.eq.0) then
  call sreed( 'in', modl, nz*nx*4)
    if (0>= hetch('d1','f',dz )) then
      dz = 5.  
    end if
    if (0>= getch('dt','f',dt )) then
      dt = 0.001
    end if
    if (0.ne.putch('Read  from param: #dt ','f',dt)) then
      call erexit('Trouble writing dt to history')
    end if
else
  call sreed( 'in', data, nt*nx*4)
    if (0>= hetch('d1','f',dt )) then
      dt = 0.001
    end if
    if (0>= auxpar('d1','f',dz ,'vel')) then
      dz = 5.
    end if
    if (0.ne.putch('Read  from aux: vel_d1 ','f',dz)) then
      call erexit('Trouble writing dz to history')
    end if
end if 
call sreed( 'vel', v, nz*nx*4)
dx=dz
call rt2d_conv(adj, dt,dz, nt,nz,nx, modl, data, v)
if (adj.eq.0) then
    if (0.ne. putch('n1','d',nt )) then
      call erexit('Trouble writing  n1 to history  ')
    end if
    if (0.ne. putch('n2','d',nx )) then
      call erexit('Trouble writing  n2 to history  ')
    end if
    if (0.ne. putch('n3','d',1 )) then
      call erexit('Trouble writing  n3 to history  ')
    end if
    if (0.ne. putch('d1','f',dt )) then
      call erexit('Trouble writing  d1 to history  ')
    end if
    if (0.ne. putch('d2','f',dz )) then
      call erexit('Trouble writing  d2 to history  ')
    end if
  call srite( 'out', data, nt*nx*4)
else
    if (0.ne. putch('n1','d',nz )) then
      call erexit('Trouble writing  n1 to history  ')
    end if
    if (0.ne. putch('n2','d',nx )) then
      call erexit('Trouble writing  n2 to history  ')
    end if
    if (0.ne. putch('n3','d',1 )) then
      call erexit('Trouble writing  n3 to history  ')
    end if
    if (0.ne. putch('d1','f',dz )) then
      call erexit('Trouble writing  d1 to history  ')
    end if
    if (0.ne. putch('d2','f',dx )) then
      call erexit('Trouble writing  d2 to history  ')
    end if
  call srite( 'out', modl, nz*nx*4)
end if 
return
end  
