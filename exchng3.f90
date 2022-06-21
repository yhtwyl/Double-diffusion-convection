subroutine exchng3(a_, sx, ex, sy, ey, sz, ez, comm3d, stridetype_EW, stridetype_NS, &
     stridetype_TB, nbrwest, nbreast, nbrsouth, nbrnorth, nbrbottom,nbrtop)
  use mpi
  
  integer sx, ex, sy, ey, sz, ez, stridetype_EW, stridetype_NS, stridetype_TB
  real a_(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1)
  integer nbrwest, nbreast, nbrsouth, nbrnorth, nbrtop, nbrbottom, comm3d
  integer ierr, ntb_,new_,nns_
!  character*4 who
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ntb_ = (ex-sx+3)*(ey-sy+3)
  !new_ = (ey-sy+1)*(ez-sz+1)
  new_ = 1
  nns_ = 1
  !nns_ = (ez-sz+1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
!  print *, nbrwest,nbreast,nbrsouth,nbrnorth,nbrbottom,nbrtop
!   if (who .eq. 'Temp') then
! !     print *, nx,nbrleft,nbrright,nbrbottom,nbrtop
!   !  These are just like the 1-d versions, except for less data
!      print *, who
!      do j = sy-1,ey+1
!         do i = sx-1,ex+1
!            print *,i,j,a_(i,j)
!         enddo
!      enddo
  !   endif
  !! For exchanging data along Top and Bottom directions
  call MPI_SENDRECV(a_(sx-1,sy-1,ez),  ntb_, MPI_real,    &
       nbrtop, 0,                              &
       a_(sx-1,sy-1,sz-1), ntb_, MPI_real,   &
       nbrbottom, 0, comm3d, MPI_STATUS_IGNORE, ierr)
  call MPI_SENDRECV(a_(sx-1,sy-1,sz),  ntb_, MPI_real,   &
       nbrbottom, 1,                          &
       a_(sx-1,sy-1,ez+1), ntb_, mpi_real,  &
       nbrtop, 1, comm3d, MPI_STATUS_IGNORE, ierr)
  !
  ! For exchanging data along east and west direction
!  print *,'exchange',sx,ex,sy,ey,sz,ez
  call MPI_SENDRECV(a_(ex,sy-1,sz-1),  new_, stridetype_EW, nbreast, 0, &
       a_(sx-1,sy-1,sz-1), new_, stridetype_EW, nbrwest, 0, &
       comm3d, MPI_STATUS_IGNORE, ierr)
  
  call MPI_SENDRECV(a_(sx,sy-1,sz-1),  new_, stridetype_EW, nbrwest, 1,&
       a_(ex+1,sy-1,sz-1), new_, stridetype_EW, nbreast, 1,&
       comm3d, MPI_STATUS_IGNORE, ierr)

!  For exchanging data along north and south directions
  call MPI_SENDRECV(a_(sx-1,ey,sz-1), nns_, stridetype_NS, nbrnorth, 0, &
       a_(sx-1,sy-1,sz-1), nns_, stridetype_NS, nbrsouth, 0, &
       comm3d, MPI_STATUS_IGNORE, ierr)
  call MPI_SENDRECV(a_(sx-1,sy,sz-1), new_, stridetype_NS, nbrsouth, 1, &
       a_(sx-1,ey+1,sz-1), nns_, stridetype_NS, nbrnorth, 1, &
       comm3d, MPI_STATUS_IGNORE, ierr)

  return
  end subroutine exchng3
