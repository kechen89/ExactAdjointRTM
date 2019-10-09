! -------------------------------------------------------------------
!
! Copyright (c) 2009 by the Society of Exploration Geophysicists.
! For more information, go to http://software.seg.org/2009/0003 .
! You must read and accept usage terms at:
! http://software.seg.org/disclaimer.txt before use.
!
! -------------------------------------------------------------------
!
!     1-Dimension
!        forward-time finite-difference extrapolation modeling ( adj=0 )
!     or the conventional RTM ( adj=1 )
!
!     Input :
!           adj      : operation type indicator ( 0 for forward / 1 for adjoint )
!           modl(nz) [when adj = 0] : model vector array
!       or  data(nt) [when adj = 1] : data vector array
!           v(nz)    : velocity model
!           nt,dt    : number of time samples and time sample distance
!           nz,dz    : number of depth samples and depth sample distance
!           snap     : switch for snap-shots output
!                         ( 0 for without / 1 for with snap-shots output)
!     Output :
!           modl(nz) [when adj = 1] : model vector array
!       or  data(nt) [when adj = 0] : data vector array
!
subroutine rt1d_conv( adj, dt,dz, nt,nz, modl, data, v, snap)
integer    adj, it,iz, nt,nz, snap
real       modl(nz), data(nt), v(nz)
real       dt,dz, dtodz2, a, a1,a2,a3, bb
real, allocatable, dimension (:,:) :: tmodl
real, allocatable, dimension (:) :: bndr
allocate (tmodl(nz,3))
allocate (bndr(20))
!
! clear data() or modl() according to adj=[0/1]
!
call adjnull(adj,0,modl,nz,data,nt)
call zero(3*nz,tmodl)
!
! Absorbing Boundary coefficients
!
do iz = 1, 20
  bb = (0.015*(20-iz))**2.
  bndr(iz) = exp(-bb)
  bndr(iz) = bndr(iz)**10
end do 
dtodz2= (dt/dz)**2.
if (adj.eq.0) then
!================== forward-time modeling starts here ======
  do iz = 1, nz  
    tmodl(iz,1) = 0.
    tmodl(iz,2) = modl(iz)
  end do 
  data(1) = modl(1)
!
! time-extrapolation loop
!
  do it = 2, nt  
!
! depth axis loop
!
    do iz = 2, nz-1  
      a = v(iz)*v(iz)*dtodz2
      tmodl(iz,3) = a*(tmodl(iz+1,2)- 2.*tmodl(iz,2)+ tmodl(iz-1,2))-&
        & tmodl(iz,1)+ 2.*tmodl(iz,2)
    end do 
!
! time-extrapolation at boundary iz=1
!
    a = v(1)*v(1)*dtodz2
    tmodl(1,3) = a*(tmodl(2,2)- 2.*tmodl(1,2))- tmodl(1,1)+ 2.*tmodl(1&
      &,2)
!
! time-extrapolation at boundary iz=nz
!
    a = v(nz)*v(nz)*dtodz2
    tmodl(nz,3) = a*(tmodl(nz-1,2)- 2.*tmodl(nz,2))- tmodl(nz,1)+ 2.&
      &*tmodl(nz,2)
    data(it) = tmodl(1,3)
    if (snap.eq.1) then
      call srite('snapout',tmodl(1,3),nz*4)
    end if
!
! Absorbing Boundary condition
!
    do iz = nz-19, nz  
      tmodl(iz,3) = bndr(nz-iz+1)*tmodl(iz,3)
    end do 
    do iz = nz-19, nz  
      tmodl(iz,2) = bndr(nz-iz+1)*tmodl(iz,2)
    end do 
    do iz = 1, nz  
      tmodl(iz,1) = tmodl(iz,2)
      tmodl(iz,2) = tmodl(iz,3)
    end do
  end do 
!================== forward-time modeling ends here ======
else
!================== reverse-time modeling starts here ====
  tmodl(1,1) = data(nt)
  tmodl(1,2) = data(nt-1)
!
! time-extrapolation loop
!
  do it = nt-2, 1, -1  
    tmodl(1,3) = data(it)
!
! Absorbing Boundary condition
!
    do iz = nz-19, nz  
      tmodl(iz,1) = bndr(nz-iz+1)*tmodl(iz,1)
    end do 
    do iz = nz-19, nz  
      tmodl(iz,2) = bndr(nz-iz+1)*tmodl(iz,2)
    end do 
!
! depth axis loop
!
    do iz = 2, nz-1  
      a = v(iz)*v(iz)*dtodz2
      tmodl(iz,3) = a*(tmodl(iz+1,2)- 2.*tmodl(iz,2)+ tmodl(iz-1,2))-&
        & tmodl(iz,1)+ 2.*tmodl(iz,2)
    end do 
!
! time extrapolation at boundary iz=nz
!
    a = v(nz)*v(nz)*dtodz2
    tmodl(nz,3) = a*(tmodl(nz-1,2)- 2.*tmodl(nz,2))- tmodl(nz,1)+ 2.&
      &*tmodl(nz,2)
    if (snap.eq.1) then
      call srite('snapout',tmodl(1,3),nz*4)
    end if
    do iz = 1, nz  
      tmodl(iz,1) = tmodl(iz,2)
      tmodl(iz,2) = tmodl(iz,3) 
      tmodl(iz,3) = 0.
    end do
  end do 
  do iz = 1, nz  
    modl(iz) = tmodl(iz,2)
  end do 
!================== reverse-time modeling ends here ====
end if 
return
end  
