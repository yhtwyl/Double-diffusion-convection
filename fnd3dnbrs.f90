!
! This routine show how to determine the neighbors in a 2-d decomposition of
! the domain. This assumes that MPI_Cart_create has already been called
!
subroutine fnd3dnbrs( comm3d, nbreast, nbrwest, nbrnorth, nbrsouth, nbrtop, nbrbottom )
  integer comm3d, nbreast, nbrwest, nbrnorth, nbrsouth, nbrtop, nbrbottom
  !     
  integer ierr,myid_
  !
  call MPI_Cart_shift( comm3d, 0,  1, nbrwest,   nbreast, ierr )
  call MPI_Cart_shift( comm3d, 1,  1, nbrsouth, nbrnorth,   ierr )
  call MPI_Cart_shift( comm3d, 2,  1, nbrbottom, nbrtop,   ierr )
  call MPI_COMM_RANK(comm3d, myid_, ierr)
  !
  !print *, nbrwest,nbreast,nbrsouth,nbrnorth,nbrbottom,nbrtop,myid_
  return
end subroutine fnd3dnbrs
!     
subroutine fnd3ddecomp( comm3d, my_id__,nx, ny, nz, sx, ex, sy, ey, sz, ez)
  integer comm3d
  integer nx, ny, nz, sx, ex, sy, ey, sz, ez
  integer dims(3), coords(3), ierr, my_id__
  logical periods(3)
  !
  call MPI_Cart_get( comm3d, 3, dims, periods, coords, ierr )
  call MPE_DECOMP1D( nx, dims(1), coords(1), sx, ex )
  call MPE_DECOMP1D( ny, dims(2), coords(2), sy, ey )
  call MPE_DECOMP1D( nz, dims(3), coords(3), sz, ez )
  !print *, 'get the diimensions', 
  return
end subroutine fnd3ddecomp
