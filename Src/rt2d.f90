! -------------------------------------------------------------------
!
! Copyright (c) 2009 by the Society of Exploration Geophysicists.
! For more information, go to http://software.seg.org/2009/0003 .
! You must read and accept usage terms at:
! http://software.seg.org/disclaimer.txt before use.
!
! -------------------------------------------------------------------
!
!     2-Dimension
!        forward-time finite-difference extrapolation modeling ( adj=0 )
!     or RTM using exact adjoint ( adj=1 )
!
!     Input :
!           adj      : operation type indicator ( 0 for forward / 1 for adjoint )
!           modl(nz,nx) [when adj = 0] : model vector array
!       or  data(nt,nx) [when adj = 1] : data vector array
!           v(nz)    : velocity model
!           nt,dt    : number of time samples and time sample distance
!           nz,dz    : number of depth samples and depth sample distance
!           nx,dx    : number of surface samples and surface sample distance
!           snap     : switch for snap-shots output
!                         (0 for without / 1 for with snap-shots output)
!
!     Output :
!           modl(nz,nx) [when adj = 1] : model vector array
!       or  data(nt,nx) [when adj = 0] : data vector array
!
subroutine rt2d( adj, dt,dx, nt,nz,nx, modl, data, v)
integer     adj, nt,nz,nx, it,iz,ix,is, bz,cz,ez, bbz
real   modl(nz,nx), data(nt,nx), v(nz,nx)
real   a,a1,a2,a3,a4, dt,dx, dtodx, bb
real, allocatable, dimension (:,:,:) :: tmodl
real, allocatable, dimension (:,:) :: v3
real, allocatable, dimension (:) :: bndr
allocate (tmodl(nz,nx,3))
allocate (v3(nz,nx))
allocate (bndr(20))
!
! clear modl() or data() according to adj=[0/1]
!
call adjnull(adj,0,modl,nz*nx, data,nt*nx)
call zero(3*nz*nx, tmodl)
!
! pre-compute damping boundary coefficients
!
do iz = 1, 20  
  bb = (0.015*(20-iz))**2.
  bndr(iz) = exp(-bb)
  bndr(iz) = bndr(iz)**10
end do 
!
! pre-compute parameters which used repeatedly
!
dtodx = dt*dt/(dx*dx)
do ix = 1, nx
  do iz = 1, nz
    v3(iz,ix) = v(iz,ix)*v(iz,ix)*dtodx
  end do 
!
end do 
if (adj.eq.0) then
!=============== forward-time modeling starts here ==========
!
!  ==> Operator [I 0 0 0 ...0 ]' in Equation (7)
!
  do iz = 1, nz
    do ix = 1, nx  
      tmodl(iz,ix,2) = modl(iz,ix)
    end do
  end do 
  do ix = 1, nx  
    data( 1,ix)   = modl( 1,ix)
  end do 
!
! time extrapolaton loop
!
  do it = 2, nt  
!
!   time extrapolaton for the range iz: 2 ~ nz-1
!   ==> Operator [0...0 -I T 0...0] in Equation (7)
!
    do iz = 2, nz-1  
      a = v3(iz,1)
      tmodl(iz,1,3) = (2.-4.*a)*tmodl(iz,1,2) - tmodl(iz,1,1) + a*(tmodl&
        &(iz,2,2) + tmodl(iz+1,1,2) + tmodl(iz-1,1,2))
      do ix = 2, nx-1  
        a = v3(iz,ix)
        tmodl(iz,ix,3) =(2.-4.*a)*tmodl(iz,ix,2) - tmodl(iz,ix,1) + a*&
          &(tmodl(iz,ix+1,2) + tmodl(iz,ix-1,2) + tmodl(iz+1,ix,2) +&
          & tmodl(iz-1,ix,2))
      end do 
      a = v3(iz,nx)
      tmodl(iz,nx,3) =(2.-4.*a)*tmodl(iz,nx,2) - tmodl(iz,nx,1) + a*&
        &(tmodl(iz,nx-1,2) + tmodl(iz+1,nx,2) + tmodl(iz-1,nx,2))
    end do 
!
! time extrapolaton at boundaries iz = 1 and iz = nx
!
    do ix = 2, nx-1  
      a = v3(1,ix)
      tmodl(1,ix,3) = (2.-4.*a)*tmodl(1,ix,2) -tmodl(1,ix,1) + a*(tmodl&
        &(2,ix,2) +tmodl(1,ix+1,2) +tmodl(1,ix-1,2))
      a = v3(nz,ix)
      tmodl(nz,ix,3) =(2.-4.*a)*tmodl(nz,ix,2) -tmodl(nz,ix,1) + a*(tmodl&
        &(nz-1,ix,2) +tmodl(nz,ix+1,2) +tmodl(nz,ix-1,2))
    end do 
!
! time extrapolaton at four corners (iz,ix) : (1,1) (nz,1) (1,nx) (nz,nx)
!
    a = v3(1,1)
    tmodl(1,1,3) =  (2.-4.*a)*tmodl(1,1,2) -tmodl(1,1,1) + a*(tmodl(2,1&
      &,2) +tmodl(1,2,2))
    a = v3(nz,1)
    tmodl(nz,1,3) = (2.-4.*a)*tmodl(nz,1,2) -tmodl(nz,1,1) + a*(tmodl&
      &(nz,2,2) +tmodl(nz-1,1,2))
    a = v3(1,nx)
    tmodl(1,nx,3) = (2.-4.*a)*tmodl(1,nx,2) -tmodl(1,nx,1) + a*(tmodl(1&
      &,nx-1,2) +tmodl(2,nx,2))
    a = v3(nz,nx)
    tmodl(nz,nx,3) =(2.-4.*a)*tmodl(nz,nx,2) -tmodl(nz,nx,1) + a*(tmodl&
      &(nz-1,nx,2) +tmodl(nz,nx-1,2))
!
! write the i-th time extrapolated wavefield onto output array data(it,ix)
!  ==> Operator S_i in Equation (7) and (8)
!
    do iz = 1, nz
      do ix = 1, nx  
        tmodl(iz,ix,1) = tmodl(iz,ix,2)
        tmodl(iz,ix,2) = tmodl(iz,ix,3)
      end do 
!
!   Applying absorbing Boundary condition to three numerical sides
!   except surface
!
    end do 
    do ix = 1, 20
      do iz = 1, nz-20  
        tmodl(iz,ix,1) = bndr(ix)*tmodl(iz,ix,1)
        tmodl(iz,ix,2) = bndr(ix)*tmodl(iz,ix,2)
      end do
    end do 
    do ix = nx-19, nx
      do iz = 1, nz-20  
        tmodl(iz,ix,1) = bndr(nx-ix+1)*tmodl(iz,ix,1)
        tmodl(iz,ix,2) = bndr(nx-ix+1)*tmodl(iz,ix,2)
      end do
    end do 
    do ix = 1, nx
      do iz = nz-19, nz  
        tmodl(iz,ix,1) = bndr(nz-iz+1)*tmodl(iz,ix,1)
        tmodl(iz,ix,2) = bndr(nz-iz+1)*tmodl(iz,ix,2)
      end do
    end do 
    do ix = 1, nx
      data(it,ix) = tmodl(1,ix,2)
    end do
  end do 
!
!   Applying absorbing Boundary condition to three numerical sides
!   except surface
!
  do it = 1, nt
    do ix = 1, 20
      data(it,ix) = bndr(ix)*data(it,ix)
    end do
  end do 
  do it = 1, nt
    do ix = nx-19, nx
      data(it,ix) = bndr(nx-ix+1)*data(it,ix)
    end do 
!=============== forward-time modeling ends here ==========
  end do
else
!=============== reverse-time migration starts here =======
  do ix = 1, nx  
    tmodl(1,ix,1) = data(nt,  ix)
    tmodl(1,ix,2) = data(nt-1,ix)
    tmodl(1,ix,3) = data(nt-2,ix)
  end do 
  cz = 3
!
! time extrapolaton loop
!
  do it = nt-1, 1, -1  
!
    cz = cz+1
    bz = min(cz,nz)
    bbz = min(bz,nz-20)
!
!   Applying absorbing Boundary condition to three numerical sides
!   except surface
!
    do ix = 1, 20
      do iz = 1, bbz  
        tmodl(iz,ix,1) = bndr(ix)*tmodl(iz,ix,1)
        tmodl(iz,ix,2) = bndr(ix)*tmodl(iz,ix,2)
      end do
    end do 
    do ix = nx-19, nx
      do iz = 1, bbz  
        tmodl(iz,ix,1) = bndr(nx-ix+1)*tmodl(iz,ix,1)
        tmodl(iz,ix,2) = bndr(nx-ix+1)*tmodl(iz,ix,2)
      end do
    end do 
    if ( bbz >= nz-19 ) then
      do ix = 1, nx
        do iz = nz-19, bbz  
          tmodl(iz,ix,1) = bndr(nz-iz+1)*tmodl(iz,ix,1)
          tmodl(iz,ix,2) = bndr(nz-iz+1)*tmodl(iz,ix,2)
        end do
      end do
    end if
!
! compute meaningful grid depth which is affected by any
! surface boundary wavefield so as to save computing time
! by avoiding unnecessary computation
!
    if ( bz .eq. nz ) then
      ez = nz - 1
    else
      ez = bz
    end if 
! time extrapolaton between iz = 2 to iz = bz
!
!   ==> Operator [ 0...0 I 0 -I 0...0] in Equation (9)
!
    do iz = 1, bz
      do ix = 1, nx
        tmodl(iz,ix,3) = tmodl(iz,ix,3) - tmodl(iz,ix,1)
      end do 
!
!   ==> Operator [ 0...0 I T' 0...0] in Equation (9)
!
    end do 
    do iz = 2, ez
      do ix = 2, nx-1  
        a  = v3(iz,  ix)
        a1 = v3(iz-1,ix)
        a2 = v3(iz+1,ix)
        a3 = v3(iz,  ix-1)
        a4 = v3(iz,  ix+1)
        tmodl(iz,ix,2) = (2.-4.*a)*tmodl(iz,ix,1) + tmodl(iz,ix,2) + a4&
          &*tmodl(iz,ix+1,1) + a3*tmodl(iz,ix-1,1) + a2*tmodl(iz+1,ix&
          &,1) + a1*tmodl(iz-1,ix,1)
      end do 
!
! time extrapolaton at boundary iz = 1
!
    end do 
    do ix = 2, nx-1  
      a  = v3(1,ix)
      a2 = v3(2,ix)
      a3 = v3(1,ix-1)
      a4 = v3(1,ix+1)
      tmodl(1,ix,2) = (2.-4.*a)*tmodl(1,ix,1) + tmodl(1,ix,2) + a4*tmodl&
        &(1,ix+1,1) + a3*tmodl(1,ix-1,1) + a2*tmodl(2,ix,1)
!
! time extrapolaton at boundary iz = nz
!
      if ( bz .eq. nz ) then
        a  = v3(nz,  ix)
        a1 = v3(nz-1,ix)
        a3 = v3(nz,  ix-1)
        a4 = v3(nz,  ix+1)
        tmodl(nz,ix,2) = (2.-4.*a)*tmodl(nz,ix,1) + tmodl(nz,ix,2) + a4&
          &*tmodl(nz,ix+1,1) + a3*tmodl(nz,ix-1,1) + a1*tmodl(nz-1,ix,1)
      end if
    end do 
!
!  extrapolaton at corner (iz,ix)=(1,1)
!
    a  = v3(1,1)
    a2 = v3(2,1)
    a4 = v3(1,2)
    tmodl(1,1,2) = (2.-4.*a)*tmodl(1,1,1) + tmodl(1,1,2) + a4*tmodl(1,2&
      &,1) + a2*tmodl(2,1,1)
!
!  extrapolaton at corner (iz,ix)=(nz,1)
!
    if ( bz .eq. nz ) then
      a  = v3(nz,  1)
      a1 = v3(nz-1,1)
      a4 = v3(nz,  2)
      tmodl(nz,1,2) = (2.-4.*a)*tmodl(nz,1,1) + tmodl(nz,1,2) + a4*tmodl&
        &(nz,2,1) + a1*tmodl(nz-1,1,1)
    end if
    do iz = 2, ez 
      a  = v3(iz,  1)
      a1 = v3(iz-1,1)
      a2 = v3(iz+1,1)
      a4 = v3(iz,2)
      tmodl(iz,1,2) = (2.-4.*a)*tmodl(iz,1,1) + tmodl(iz,1,2) + a4*tmodl&
        &(iz,2,1) + a2*tmodl(iz+1,1,1) + a1*tmodl(iz-1,1,1)
      a  = v3(iz,  nx)
      a1 = v3(iz-1,nx)
      a2 = v3(iz+1,nx)
      a3 = v3(iz,  nx-1)
      tmodl(iz,nx,2) = (2.-4.*a)*tmodl(iz,nx,1) + tmodl(iz,nx,2) + a3&
        &*tmodl(iz,nx-1,1) + a2*tmodl(iz+1,nx,1) + a1*tmodl(iz-1,nx,1)
    end do 
!
!  extrapolaton at corner (iz,ix)=(1,nx)
!
    a  = v3(1,  nx)
    a2 = v3(1+1,nx)
    a3 = v3(1,  nx-1)
    tmodl(1,nx,2) = (2.-4.*a)*tmodl(1,nx,1) + tmodl(1,nx,2) + a3*tmodl&
      &(1,nx-1,1) + a2*tmodl(2,nx,1) 
!
!  extrapolaton at corner (iz,ix)=(nz,nx)
!
    if ( bz .eq. nz ) then
      a  = v3(nz,  nx)
      a1 = v3(nz-1,nx)
      a3 = v3(nz,  nx-1)
      tmodl(nz,nx,2) = (2.-4.*a)*tmodl(nz,nx,1) + tmodl(nz,nx,2) + a3&
        &*tmodl(nz,nx-1,1) + a1*tmodl(nz-1,nx,1)
    end if
!
! shift wavefield in temporary array for next time-extrapolation
!
    do iz = 1, nz
      do ix = 1, nx  
        tmodl(iz,ix,1) = tmodl(iz,ix,2)
        tmodl(iz,ix,2) = tmodl(iz,ix,3)
      end do 
!
! Insert the input surface boundary wavefield
!   ==> Operator S'_i in Equation (9) and (8)
!
    end do 
    if ( it > 2 ) then
      do iz = 2, nz
        do ix = 1, nx  
          tmodl(iz,ix,3) = 0.
        end do
      end do 
      do ix = 1, nx  
        tmodl( 1,ix,3) = data(it-2,ix)
      end do
    end if
  end do 
!
! write the final wavefield onto output array modl(iz,ix)
!  ==> Operator [I 0 0 0 ...0 ] in Equation (9)
!
  do iz = 1, nz
    do ix = 1, nx  
      modl(iz,ix) = tmodl(iz,ix,1)
    end do 
!=============== reverse-time migration ends here =======
  end do
end if 
return
end  
