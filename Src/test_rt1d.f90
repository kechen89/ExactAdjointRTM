! -------------------------------------------------------------------
!
! Copyright (c) 2009 by the Society of Exploration Geophysicists.
! For more information, go to http://software.seg.org/2009/0003 .
! You must read and accept usage terms at:
! http://software.seg.org/disclaimer.txt before use.
!
! -------------------------------------------------------------------
!
!  1-D Time-extrapolation modeling and Reverse-time Migration (RTM)
!  using subroutine rt1d() where RTM with accurate adjoint implementation
!
!  Input
!       adj=[0/1] : operation type indicator ( 0 for forward / 1 for adjoint )
!       in.H : SEPlib header file for 1D wavefield
!              should contains n1:nt, d1:dt
!       vel=vel.H : SEPlib header file for the 1D velocity function
!                   should contains n1:nz, d1:dz
!       snap=[0/1] : switch for snap-shots output
!              (0 for without/1 for with snap-shots output)
!       snapout=snap.H : snap-shots file output
!
!  Output
!       out.H : data [when adj=0] or image [when adj=1]
!  Usage
!
!	test_rt1d.x < in.H adj=[0/1] vel=vel.H snap=[0/1] snapout=snap.H > out.H
!
!%----------------
program zyxabc
implicit none
real, allocatable, dimension (:) :: modl
real, allocatable, dimension (:) :: v
real, allocatable, dimension (:) :: data
integer n1,n3, nt,nz, adj, snap
integer getch,hetch,putch,auxpar
call initpar()
call doc('./../Src/test_rt1d.f90')
!
! read parameters adj[=0] and snap[=0] from command line
!
  if (0>= getch('adj','d',adj )) then
    adj = 0
  end if
  if (0.ne.putch('Read  from param: #adj ','d',adj)) then
    call erexit('Trouble writing adj to history')
  end if
  if (0>= getch('snap','d',snap )) then
    snap = 0
  end if
  if (0.ne.putch('Read  from param: #snap ','d',snap)) then
    call erexit('Trouble writing snap to history')
  end if
!
! read parameters n1[=nz/nt] and nt/nz from history file
!
if (adj.eq.0) then
    if (0>= hetch('n1','d',nz )) then
      call erexit('Trouble reading  n1 from history  ')
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
    if (0>= auxpar('n1','d',nz ,'vel')) then
      call erexit('Trouble reading  n1 from aux vel ')
    end if
    if (0.ne.putch('Read  from aux: vel_n1 ','d',nz)) then
      call erexit('Trouble writing nz to history')
    end if
end if 
allocate (modl(nz))
allocate (v(nz))
allocate (data(nt))
call submain ( nt,nz, modl,data, v, adj,snap)
end program zyxabc 
!
subroutine submain( nt,nz, modl,data, v, adj,snap)
integer n1,n2,n3, nt,nz, i,j, adj,snap
real    d1,d3, dt,dz, o1,o2,o3 
real    modl(nz)
real    v(nz)
real    data(nt)
integer getch,putch,auxputch
call sreed( 'vel', v, nz*4)
if (adj.eq.0) then
  call sreed( 'in', modl, nz*4)
else
  call sreed( 'in', data, nt*4)
end if
  if (0>= getch('dt','f',dt )) then
    dt = .004
  end if
  if (0.ne.putch('Read  from param: #dt ','f',dt)) then
    call erexit('Trouble writing dt to history')
  end if
  if (0>= getch('dz','f',dz )) then
    dz = 25.
  end if
  if (0.ne.putch('Read  from param: #dz ','f',dz)) then
    call erexit('Trouble writing dz to history')
  end if
if (snap.eq.1) then
    if (0.ne. auxputch('n1','d',nz ,'snapout')) then
      call erexit('Trouble writing  n1 to aux snapout ')
    end if
    if (0.ne. auxputch('n2','d',nt ,'snapout')) then
      call erexit('Trouble writing  n2 to aux snapout ')
    end if
end if
call rt1d( adj, dt,dz, nt,nz, modl, data, v, snap)
if (adj.eq.0) then
    if (0.ne. putch('n1','d',nt )) then
      call erexit('Trouble writing  n1 to history  ')
    end if
    if (0.ne. putch('n2','d',1 )) then
      call erexit('Trouble writing  n2 to history  ')
    end if
    if (0.ne. putch('n3','d',1 )) then
      call erexit('Trouble writing  n3 to history  ')
    end if
    if (0.ne. putch('d1','f',dt )) then
      call erexit('Trouble writing  d1 to history  ')
    end if
  call srite( 'out', data, nt*4)
else
    if (0.ne. putch('n1','d',nz )) then
      call erexit('Trouble writing  n1 to history  ')
    end if
    if (0.ne. putch('n2','d',1 )) then
      call erexit('Trouble writing  n2 to history  ')
    end if
    if (0.ne. putch('n3','d',1 )) then
      call erexit('Trouble writing  n3 to history  ')
    end if
    if (0.ne. putch('d1','f',dz )) then
      call erexit('Trouble writing  d1 to history  ')
    end if
  call srite( 'out', modl, nz*4)
end if 
return
end  
