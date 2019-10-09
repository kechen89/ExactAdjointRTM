! -------------------------------------------------------------------
!
! Copyright (c) 2009 by the Society of Exploration Geophysicists.
! For more information, go to http://software.seg.org/2009/0003 .
! You must read and accept usage terms at:
! http://software.seg.org/disclaimer.txt before use.
!
! -------------------------------------------------------------------
!
program zyxabc
implicit none
real, allocatable, dimension (:,:) :: data
real, allocatable, dimension (:,:) :: modl
real, allocatable, dimension (:,:) :: v
real, allocatable, dimension (:) :: res
!  2-D wavefield inversion using Finite Difference time extrapolation operator
!
!  Input
!       in.H    : SEPlib header file for 2D wavefield
!                 should contains n1:nt, n2:nx, d1:dt, d2:dx
!       niter   : number of iteration in the conjugate gradient algorithm
!       vel=vel.H : SEPlib header file for the 2D velocity function
!                      should contains n1:nz, n2:nx, d1:dz, d2:dx
!       mode=[0/1] : 0 for using accurate adjoint operator pair for inversion
!                  : 1 for using approximate adjoint operator pair for inversion
!
!  Output
!       out.H   : inverse image --> will have the size as vel.H
!       res=res.H : SEPlib header file for residuals w.r.t. iteration,
!                   contains n1:niter+1
!
!  Usange
!
!       Rt2dinv.x niter=10 vel=vel.H res=res.H < in.H > out.H
!
! ----------------
integer n1,n2, nt,nz,nx, niter
integer getch,hetch,putch,auxpar
call initpar()
call doc('/home/jun/ADJRTM/Prog2/Src/Rt2dinv.rs')
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
  if (0>= auxpar('n2','d',nx ,'vel')) then
    call erexit('Trouble reading  n2 from aux vel ')
  end if
  if (0.ne.putch('Read  from aux: vel_n2 ','d',nx)) then
    call erexit('Trouble writing nx to history')
  end if
  if (0>= getch('niter','d',niter )) then
    niter = 20
  end if
  if (0.ne.putch('Read  from param: #niter ','d',niter)) then
    call erexit('Trouble writing niter to history')
  end if
allocate (data(nt,nx))
allocate (modl(nz,nx))
allocate (v(nz,nx))
allocate (res(niter+1))
call submain (nt,nz,nx, modl,data,v,res, niter)
end program zyxabc 
subroutine submain(nt,nz,nx, modl,data,v,res, niter)
integer n1,n2,n3, nt,nz,nx, niter, mode
real    d1,d2,d3, dt,dz,dx, o1,o2
real    data(nt,nx)
real    modl(nz,nx)
real    v(nz,nx)
real    res(niter+1)
integer getch,hetch,putch,auxpar,auxputch
  if (0>= getch('mode','d',mode )) then
    mode = 0
  end if
  if (0.ne.putch('Read  from param: #mode ','d',mode)) then
    call erexit('Trouble writing mode to history')
  end if
  if (0>= getch('dt','f',dt )) then
    dt = .001
  end if
  if (0.ne.putch('Read  from param: #dt ','f',dt)) then
    call erexit('Trouble writing dt to history')
  end if
  if (0>= getch('dz','f',dz )) then
    dz = 5.
  end if
  if (0.ne.putch('Read  from param: #dz ','f',dz)) then
    call erexit('Trouble writing dz to history')
  end if
  if (0>= hetch('d1','f',dt )) then
    dt = .001
  end if
  if (0>= hetch('d2','f',dx )) then
    dx = dz
  end if
call sreed(  'in', data, nt*nx*4)
  if (0>= auxpar('d1','f',dz ,'vel')) then
    dz = 5.
  end if
  if (0.ne.putch('Read  from aux: vel_d1 ','f',dz)) then
    call erexit('Trouble writing dz to history')
  end if
call sreed( 'vel',    v, nz*nx*4)
write(0,*) 'iter=', niter
write(0,*) 'mode=', mode
write(0,*) 'nt=', nt, 'nz=', nz, 'nx=', nx
if (mode.eq.0) then
  call rt2dinv(nz,nx,modl, nt,data, dt,dz,v, res, niter)
else
  call rt2dinv2(nz,nx,modl, nt,data, dt,dz,v, res, niter)
end if
  if (0.ne. auxputch('n1','d',niter+1 ,'res')) then
    call erexit('Trouble writing  n1 to aux res ')
  end if
call srite( 'res', res, (niter+1)*4)
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
  if (0.ne. putch('d2','f',dz )) then
    call erexit('Trouble writing  d2 to history  ')
  end if
  if (0.ne. putch('d3','f',0. )) then
    call erexit('Trouble writing  d3 to history  ')
  end if
  if (0.ne. putch('o1','f',0. )) then
    call erexit('Trouble writing  o1 to history  ')
  end if
  if (0.ne. putch('o2','f',0. )) then
    call erexit('Trouble writing  o2 to history  ')
  end if
  if (0.ne. putch('o3','f',0. )) then
    call erexit('Trouble writing  o3 to history  ')
  end if
call srite( 'out', modl, nz*nx*4)
return
end  
