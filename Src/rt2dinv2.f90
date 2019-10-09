! -------------------------------------------------------------------
!
! Copyright (c) 2009 by the Society of Exploration Geophysicists.
! For more information, go to http://software.seg.org/2009/0003 .
! You must read and accept usage terms at:
! http://software.seg.org/disclaimer.txt before use.
!
! -------------------------------------------------------------------
!
!     2-Dimensional wavefield inversion
!     for the forward-time finite-difference extrapolation modeling
!
!     Input
!           data(nt,nx) : data vector array
!           v(nz,nx)    : velocity model
!           nt,dt    : number of time samples and time sample distance
!           nz,dz    : number of depth samples and depth sample distance
!           nx,dx    : number of surface samples and surface sample distance
!           niter    : number of iteration in the conjugate gradient step
!     Output
!           modl(nz,nx) : model vector array
!           res(niter+1) : residuals w.r.t. iteration in mean square error
!
subroutine rt2dinv2( nz,nx,model, nt,data, dt,dz,vel,res, niter)
integer   nz,nx,nt, iz,ix,it,iter, niter
real      model(nz,nx),data(nt,nx), dt,dz,vel(nz,nx),res(niter+1), dot
real  dmodel(nz,nx), smodel(nz,nx), dr(nt,nx), sr(nt,nx), rr(nt,nx)
do iz= 1, nz
  do ix= 1, nx
    model(iz,ix) = 0.0
  end do 
!
end do 
do it= 1, nt
  do ix= 1, nx
    rr(it,ix) = -data(it,ix)
  end do 
!
end do 
do iter = 0, niter  
  call rt2d_conv( 1, dt,dz, nt,nz,nx, dmodel, rr, vel)
  call rt2d_conv( 0, dt,dz, nt,nz,nx, dmodel, dr, vel)
  call cgstep(iter, nz*nx, model, dmodel, smodel, nt*nx, rr, dr, sr)
  res(iter+1) = dot(nt*nx,rr,rr)/(nt*nx)
  write(0,10) 'iter=',iter,'  of', niter
  10 format(A,I3,A,I3)
end do 
return
end  
