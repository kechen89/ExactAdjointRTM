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
!     or RTM using the conventional way ( adj=1 )
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
!                        (0 for without / 1 for with snap-shots output)
!     Output :
!           modl(nz,nx) [when adj = 1] : model vector array
!       or  data(nt,nx) [when adj = 0] : data vector array
!
subroutine rt2d_conv( adj, dt,dx, nt,nz,nx, modl, data, v)
integer adj, nt,nz,nx, it,iz,ix,is, cz,bz,ez, bbz,bbbz
real    modl(nz,nx), data(nt,nx), v(nz,nx)
real    a,a1,a2,a3,a4, dt,dx, dtodx, bb
real, allocatable, dimension (:,:,:) :: tmodl
real, allocatable, dimension (:,:) :: v2
real, allocatable, dimension (:) :: bndr
allocate (tmodl(nz,nx,3))
allocate (v2(nz,nx))
allocate (bndr(20))
!
! clear modl() or data() according to adj=[0/1]
!
call adjnull(adj,0,modl,nz*nx,data,nt*nx)
call zero(3*nx*nz,tmodl)
!
! damping boundary coefficients
!
do iz = 1, 20
  bb = (0.015*(20-iz))**2.
  bndr(iz) = exp(-bb)
  bndr(iz) = bndr(iz)**10
end do 
dtodx = dt/dx
do ix = 1, nx
  do iz = 1, nz
    v2(iz,ix) = v(iz,ix)*v(iz,ix)*dtodx*dtodx
  end do
end do 
if (adj.eq.0) then
!=============== forward-time modeling starts here ==========
  do iz = 1, nz
    do ix = 1, nx  
      tmodl(iz,ix,1) = 0.
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
! depth axis from iz=2 to iz=nz-1
!
    do iz = 2, nz-1  
      a = v2(iz,1)
      tmodl(iz,1,3) = (2.-4.*a)*tmodl(iz,1,2)- tmodl(iz,1,1)+ a *(tmodl&
        &(iz,2,2)+ tmodl(iz+1,1,2)+ tmodl(iz-1,1,2))
      do ix = 2, nx-1  
        a = v2(iz,ix)
        tmodl(iz,ix,3) =(2.-4.*a)*tmodl(iz,ix,2)- tmodl(iz,ix,1)+ a *&
          &(tmodl(iz,ix+1,2)+ tmodl(iz,ix-1,2)+ tmodl(iz+1,ix,2)+ tmodl&
          &(iz-1,ix,2))
      end do 
      a = v2(iz,nx)
      tmodl(iz,nx,3) =(2.-4.*a)*tmodl(iz,nx,2)- tmodl(iz,nx,1)+ a *(tmodl&
        &(iz,nx-1,2)+ tmodl(iz+1,nx,2)+ tmodl(iz-1,nx,2))
    end do 
!
! time extrapolaton at boundaries iz=1 and iz=nx
!
    do ix = 2, nx-1  
      a = v2(1,ix)
      tmodl(1,ix,3) = (2.-4.*a)*tmodl(1,ix,2)- tmodl(1,ix,1)+ a *(tmodl&
        &(2,ix,2)+ tmodl(1,ix+1,2)+ tmodl(1,ix-1,2))
      a = v2(nz,ix)
      tmodl(nz,ix,3) =(2.-4.*a)*tmodl(nz,ix,2)- tmodl(nz,ix,1)+ a *(tmodl&
        &(nz-1,ix,2)+ tmodl(nz,ix+1,2)+ tmodl(nz,ix-1,2))
    end do 
!
! time extrapolaton at four corners (iz,ix):(1,1)(nz,1)(1,nx)(nz,nx)
!
    a = v2(1,1)
    tmodl(1,1,3) =  (2.-4.*a)*tmodl(1,1,2)- tmodl(1,1,1)+ a *(tmodl(2,1&
      &,2)+ tmodl(1,2,2))
    a = v2(nz,1)
    tmodl(nz,1,3) = (2.-4.*a)*tmodl(nz,1,2)- tmodl(nz,1,1)+ a *(tmodl&
      &(nz,2,2)+ tmodl(nz-1,1,2))
    a = v2(1,nx)
    tmodl(1,nx,3) = (2.-4.*a)*tmodl(1,nx,2)- tmodl(1,nx,1)+ a *(tmodl(1&
      &,nx-1,2)+ tmodl(2,nx,2))
    a = v2(nz,nx)*dtodx
    tmodl(nz,nx,3) =(2.-4.*a)*tmodl(nz,nx,2)- tmodl(nz,nx,1)+ a *(tmodl&
      &(nz-1,nx,2)+ tmodl(nz,nx-1,2))
!
! Absorbing Boundary condition
!
    do ix = 1, 20
      do iz = 1, nz-20  
        tmodl(iz,ix,3) = bndr(ix)*tmodl(iz,ix,3)
        tmodl(iz,ix,2) = bndr(ix)*tmodl(iz,ix,2)
      end do
    end do 
    do ix = nx-19, nx
      do iz = 1, nz-20  
        tmodl(iz,ix,3) = bndr(nx-ix+1)*tmodl(iz,ix,3)
        tmodl(iz,ix,2) = bndr(nx-ix+1)*tmodl(iz,ix,2)
      end do
    end do 
    do ix = 1, nx
      do iz = nz-19, nz  
        tmodl(iz,ix,3) = bndr(nz-iz+1)*tmodl(iz,ix,3)
        tmodl(iz,ix,2) = bndr(nz-iz+1)*tmodl(iz,ix,2)
      end do 
!
! write the i-th extrapolated wavefield onto output array data(it,ix)
!
    end do 
    do ix = 1, nx
      data(it,ix) = tmodl(1,ix,3)
    end do 
    do iz = 1, nz
      do ix = 1, nx  
        tmodl(iz,ix,1) = tmodl(iz,ix,2)
        tmodl(iz,ix,2) = tmodl(iz,ix,3)
      end do
    end do
  end do 
!=============== forward-time modeling ends here ==========
else
!=============== reverse-time migration starts here =======
  do ix = 1, nx  
    tmodl(1,ix,1) = data(nt,  ix)
    tmodl(1,ix,2) = data(nt-1,ix)
  end do 
  cz = 3
!
! time extrapolaton loop
!
  do it = nt-2, 1, -1  
! compute meaningful grid depth which is affected by any
! surface boundary wavefield so as to save computing time
! by avoiding unnecessary computation
    cz = cz + 1
    bz = min(cz,nz)
    if ( bz .eq. nz ) then
      ez = nz - 1
    else
      ez = bz
    end if 
    bbz = min(bz,nz-20)
!
! Absorbing boundary condition
!
    do ix = 1, 20
      do iz = 1, bbz  
        tmodl(iz,ix,2) = bndr(ix)*tmodl(iz,ix,2)
        tmodl(iz,ix,1) = bndr(ix)*tmodl(iz,ix,1)
      end do
    end do 
    do ix = nx-19, nx
      do iz = 1, bbz  
        tmodl(iz,ix,2) = bndr(nx-ix+1)*tmodl(iz,ix,2)
        tmodl(iz,ix,1) = bndr(nx-ix+1)*tmodl(iz,ix,1)
      end do
    end do 
    if ( bbz >= nz-19) then
      bbbz = min(bbz,nz)
      do ix = 1, nx
        do iz = nz-19, bbbz  
          tmodl(iz,ix,2) = bndr(nz-iz+1)*tmodl(iz,ix,2)
          tmodl(iz,ix,1) = bndr(nz-iz+1)*tmodl(iz,ix,1)
        end do
      end do
    end if
!
! time extrapolaton from iz = 2 to iz = ez
!
    do iz = 2, ez  
      a = v2(iz,1)
      tmodl(iz,1,3) = (2.-4.*a)*tmodl(iz,1,2)- tmodl(iz,1,1)+ a *(tmodl&
        &(iz,2,2)+ tmodl(iz+1,1,2)+ tmodl(iz-1,1,2))
      do ix = 2, nx-1  
        a = v2(iz,ix)
        tmodl(iz,ix,3) =(2.-4.*a)*tmodl(iz,ix,2)- tmodl(iz,ix,1)+ a *&
          &(tmodl(iz,ix+1,2)+ tmodl(iz,ix-1,2)+ tmodl(iz+1,ix,2)+ tmodl&
          &(iz-1,ix,2))
      end do 
      a = v2(iz,nx)
      tmodl(iz,nx,3) =(2.-4.*a)*tmodl(iz,nx,2)- tmodl(iz,nx,1)+ a *(tmodl&
        &(iz,nx-1,2)+ tmodl(iz+1,nx,2)+ tmodl(iz-1,nx,2))
    end do 
!
! time extrapolaton at boundary iz=nz
!
    if ( bz.eq.nz) then
      do ix = 2, nx-1  
        a = v2(nz,ix)
        tmodl(nz,ix,3) =(2.-4.*a)*tmodl(nz,ix,2)- tmodl(nz,ix,1)+ a *&
          &(tmodl(nz-1,ix,2)+ tmodl(nz,ix+1,2)+ tmodl(nz,ix-1,2))
      end do 
      a = v2(nz,1)
      tmodl(nz,1,3) = (2.-4.*a)*tmodl(nz,1,2)- tmodl(nz,1,1)+ a *(tmodl&
        &(nz,2,2)+ tmodl(nz-1,1,2))
      a = v2(nz,nx)
      tmodl(nz,nx,3) =(2.-4.*a)*tmodl(nz,nx,2)- tmodl(nz,nx,1)+ a *(tmodl&
        &(nz-1,nx,2)+ tmodl(nz,nx-1,2))
    end if
    do ix = 1, nx  
      tmodl( 1,ix,3) = data(it,ix)
    end do 
!
! shift wavefield in temporary array for next time-extrapolation
!
    do iz = 1, nz
      do ix = 1, nx  
        tmodl(iz,ix,1) = tmodl(iz,ix,2)
        tmodl(iz,ix,2) = tmodl(iz,ix,3)
      end do
    end do
  end do 
!
! Insert the input surface boundary wavefield
!
  do iz = 1, nz
    do ix = 1, nx  
      modl(iz,ix) = tmodl(iz,ix,2)
    end do 
!=============== reverse-time migration ends here =======
  end do
end if 
return
end  
