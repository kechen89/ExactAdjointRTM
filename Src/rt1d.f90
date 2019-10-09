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
!     or reverse-time finite-difference extrapolation migration ( adj=1 )
!
!     Input :
!           adj      : operation type indicator ( 0 for forward / 1 for adjoint )
!           modl(nz) [when adj = 0] : model vector array
!       or  data(nt) [when adj = 1] : data vector array
!           v(nz)    : velocity model
!           nt,dt    : number of time samples and time sample distance
!           nz,dz    : number of depth samples and depth sample distance
!           snap     : switch for snap-shots output
!                      ( 0 for without / 1 for with snap-shots output)
!
!     Output :
!           modl(nz) [when adj = 1] : model vector array
!       or  data(nt) [when adj = 0] : data vector array
!
subroutine rt1d( adj, dt,dz, nt,nz, modl, data, v, snap)
integer          adj, nt,nz, it,iz, snap
real                  dt,dz, dtz, a, a1,a2,a3, bb
real                  modl(nz), data(nt), v(nz)
real, allocatable, dimension (:,:) :: tm
real, allocatable, dimension (:) :: bndr
allocate (tm(nz,3))
allocate (bndr(20))
!
! clear data() or modl() according to adj=[0/1]
!
call adjnull(adj, 0, modl,nz, data,nt)
!
! clear the temporary array
!
call zero(3*nz, tm)
!
! Absorbing Boundary coefficients
!
do iz = 1, 20
  bb = (0.015*(20-iz))**2.
  bndr(iz) = exp(-bb)
  bndr(iz) = bndr(iz)**10
end do 
!
! forward-time modeling
!
if (adj.eq.0) then
!============== forward-time modeling starts here ======
!
!   ==> Operator [I 0 0 0 ...0 ]' in Equation (7)
!
  do iz = 1, nz  
    tm(iz,1) = 0.
    tm(iz,2) = modl(iz)
  end do 
  data(1) = modl(1)
!
! time-extrapolaton loop
!
  do it = 2, nt  
!
! Absorbing Boundary condition
!
    do iz = nz-19, nz  
      tm(iz,1) = bndr(nz-iz+1)*tm(iz,1)
    end do 
    do iz = nz-19, nz  
      tm(iz,2) = bndr(nz-iz+1)*tm(iz,2)
    end do 
!
! depth axis from iz = 2 to iz = nz-1
!   ==> Operator [-I T] in Equation (5)
!
    do iz = 2, nz-1  
      a = (v(iz)*dt/dz)**2.
      tm(iz,3) = a*(tm(iz+1,2) - 2.*tm(iz,2)+tm(iz-1,2))- tm(iz,1) + 2.&
        &*tm(iz,2)
    end do 
!
! time extrapolaton at boundary iz = 1
!
    a = (v(1)*dt/dz)**2.
    tm(1,3) = a*(tm(2,2) - 2.*tm(1,2)) - tm(1,1) + 2.*tm(1,2)
!
! time extrapolaton at boundary iz = nz
!
    a = (v(nz)*dt/dz)**2.
    tm(nz,3) = a*(tm(nz-1,2) - 2.*tm(nz,2))- tm(nz,1) + 2.*tm(nz,2)
!
!  write the i-th extrapolated wavefield onto output array data()
!  ==> Operator S_i in Equation (8)
!
    data(it) = tm(1,3)
!
! output snapshots of extrapolated wavefield
!
    if (snap.eq.1 .and. it.ne.nt) then
      call srite('snapout',tm(1,3),nz*4)
    end if
    do iz = 1, nz  
      tm(iz,1) = tm(iz,2)
      tm(iz,2) = tm(iz,3)
    end do
  end do 
!============ forward-time modeling ends here ========
else
!============ reverse-time migration starts here =====
  dtz= (dt/dz)**2.
  tm(1,1) = data(nt)
  tm(1,2) = data(nt-1)
  tm(1,3) = data(nt-2)
!
! reverse time extrapolaton loop
!
  do it = nt-1, 1, -1  
!
! Absorbing Boundary condition
!
    do iz = nz-19, nz  
      tm(iz,1) = bndr(nz-iz+1)*tm(iz,1)
    end do 
    do iz = nz-19, nz  
      tm(iz,2) = bndr(nz-iz+1)*tm(iz,2)
    end do 
!
    do iz = 1, nz
      tm(iz,3) = tm(iz,3) - tm(iz,1)
    end do 
!
! depth axis from iz = 2 to iz = nz-1
!   ==> Operator [ 0...0 -I T 0...0 ]' in Equation (9)
!
    do iz = 2, nz-1 
      a1 = v(iz-1)*v(iz-1)*dtz
      a2 = v(iz  )*v(iz  )*dtz
      a3 = v(iz+1)*v(iz+1)*dtz
      tm(iz,2) = a1*tm(iz-1,1) - 2.*a2*tm(iz,1) + a3*tm(iz+1,1)+ 2.*tm&
        &(iz,1) + tm(iz,2)
    end do 
!
! reverse time extrapolaton at boundary iz = 1
!
    a2 = v(1)*v(1)*dtz
    a3 = v(2)*v(2)*dtz
    tm( 1,2) = a3*tm(2,1) - 2.*a2*tm(1,1) + 2.*tm(1,1) + tm(1,2)
!
! reverse time extrapolaton at boundary iz = nz
!
    a1 = v(nz-1)*v(nz-1)*dtz
    a2 = v(nz  )*v(nz  )*dtz
    tm(nz,2) = a1*tm(nz-1,1) - 2.*a2*tm(nz,1) + 2.*tm(nz,1) + tm(nz,2)
!
! output snapshots of extrapolated wavefield
!
    if (snap.eq.1 .and. it .ne. 1 ) then
      call srite('snapout',tm(1,2),nz*4)
    end if
!
! shift wavefield in temporary array for next time-extrapolation
!
    do iz = 1, nz  
      tm(iz,1) = tm(iz,2)
      tm(iz,2) = tm(iz,3)
    end do 
!
! Insert the input surface boundary wavefield
!  ==> Operator S'_i in Equation (9)
!
    if ( it > 2 ) then
      do iz = 1, nz
        tm(iz,3) = 0.
      end do 
      tm( 1,3) = data(it-2)
    end if
  end do 
!
! write the final wavefield onto output array modl()
!   ==> Operator [I 0 0 ...0 ] in Equation (9)
!
  do iz = 1, nz  
    modl(iz) = tm(iz,1)
  end do 
!=============  reverse-time migration ends here ========
end if 
return
end  
