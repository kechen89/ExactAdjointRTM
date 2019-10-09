! -------------------------------------------------------------------
!
! Copyright (c) 2009 by the Society of Exploration Geophysicists.
! For more information, go to http://software.seg.org/2009/0003 .
! You must read and accept usage terms at:
! http://software.seg.org/disclaimer.txt before use.
!
! -------------------------------------------------------------------
!
!     1-Dimensional wavefield inversion
!     for the forward-time finite-difference extrapolation modeling
!     using the conventional RTM as the adjoint operator
!
!     Input
!           data(nt) : data vector array
!           v(nz)    : velocity model
!           nt,dt    : number of time samples and time sample distance
!           nz,dz    : number of depth samples and depth sample distance
!           niter    : number of iteration in the conjugate gradient step
!     Output
!           modl(nz) : model vector array
!
subroutine rt1dinv2( nz,model, nt,data, dt,dz,vel, res, niter)
integer              nz,nt, niter,iz,it,iter
real     model(nz),data(nt), dt,dz,vel(nz), res(niter+1), dot
real  dmodel(nz), smodel(nz), dr(nt), sr(nt), rr(nt)
do iz= 1, nz
  model(iz) = 0.0
end do 
!
do it= 1, nt
  rr(it) = -data(it)
end do 
!
do iter = 0, niter  
  call rt1d_conv( 1, dt,dz, nt,nz, dmodel, rr, vel,0)
  call rt1d_conv( 0, dt,dz, nt,nz, dmodel, dr, vel,0)
  call cgstep(iter, nz, model, dmodel, smodel, nt, rr, dr, sr)
  res(iter+1) = dot(nt,rr,rr)/nt
end do 
return
end  
