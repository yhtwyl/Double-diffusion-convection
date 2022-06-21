program main
  use mpi
  integer ierr, myid,numprocs,dims(3),sx,sy,sz,ex,ey,ez,nx,ny,nz,rc,comm3d
  REAL, DIMENSION(100):: DATA_
  logical periods(3)
  double precision t1,t2
 external diff2d
  data periods/3*.false./
  
 integer(kind=MPI_OFFSET_KIND), parameter :: zero_off = 0
  
  t1 = MPI_WTIME();
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
  print *, myid,numprocs
  if (myid .eq. 0) then
     
     OPEN(UNIT=10,FILE='3_pnew.inp')
!     OPEN(UNIT=17,FILE=FILE17)
!     OPEN(UNIT=18,FILE=FILE18)
!     OPEN(UNIT=19,FILE=FILE19)
     

     !WRITE(*,*)'Supply - NX,NY,NZ,BIG1,BIG2,PHID1,PHID2,PHID3 '
     READ(10,*)
!     READ(10,*) NX,NY,BIG1,BIG2,PHID1,PHID2
     READ(10,*)  ( DATA_(I), I =1,8)
!     !WRITE(*,*) NX,NY,BIG1,BIG2,PHID1,PHID2
     !WRITE(*,*) ( DATA_(I), I =1,8)
     
     !WRITE(*,*)'- GRAST,IT,GRASM,IC,PRNDL,SCHMDT,THIGH,TLOW '
     READ(10,*)
     !     READ(10,*)   GRAST,IT,GRASM,IC,PRNDL,SCHMDT,THIGH,TLOW
     READ(10,*)  ( DATA_(8+I), I =1,8)
!     !WRITE(*,*)   GRAST,IT,GRASM,IC,PRNDL,SCHMDT,THIGH,TLOW
     !WRITE(*,*) ( DATA_(8+I), I =1,8)
     
     !WRITE(*,*)'  DTAU,TAUMAX,IFLAG,ARXZ,ARYZ,CONMAS,NACELL '
     READ(10,*)
     !     READ(10,*)   DTAU,TAUMAX,IFLAG,ARXZ,ARYZ,CONMAS,NACELL
     READ(10,*) ( DATA_(16+I), I =1,7)
!     !WRITE(*,*)   DTAU,TAUMAX,IFLAG,ARXY,CONMAS,NACELL
     !WRITE(*,*) ( DATA_(16+I), I =1,7)
     
     !WRITE(*,*)'  ITLIM,ITMASS,NTAU,NITMAX,NUSOPT,NTAUSTEP,NFUN,INIT_ITER'
     READ(10,*)  
!     READ(10,*)   ITLIM,ITMASS,NTAU,NITMAX,NUSOPT,NTAUSTEP,NFUN,NKINE, &
!          INIT_ITER
     READ(10,*) ( DATA_(23+I), I =1,9)
!     !WRITE(*,*)   ITLIM,ITMASS,NTAU,NITMAX,NUSOPT,NTAUSTEP,NFUN,NKINE, &
!          INIT_ITER
     !WRITE(*,*) ( DATA_(23+I), I =1,9)
     
     !WRITE(*,*)'  ERMASS,ERTEM,ERTEC,ERU,ERV,ERW,SIGN  '
     READ(10,*) 
!     READ(10,*)   ERMASS,ERTEM,ERTEC,ERU,ERV,ERW, SIGN
     READ(10,*) ( DATA_(32+I), I =1,7)
!     !WRITE(*,*)   ERMASS,ERTEM,ERTEC,ERU,ERV,SIGN
     !WRITE(*,*) ( DATA_(32+I), I =1,7)
     
     !WRITE(*,*)'  RELAXU, RELAXV, RELAXW, RELAXT, RELAXC  '
     READ(10,*)
     READ(10,*) ( DATA_(39+I), I =1,5)
     !WRITE(*,*) ( DATA_(39+I), I =1,5)
!     READ(10,*)  ( RELAX(I), I =1,4)
!     !WRITE(*,*)  ( RELAX(I), I =1,4)

     !WRITE(*,*)' TSTBOT,CSTBOT,TINBOT,CINBOT '
     READ(10,*)
     
     READ(10,*) ( DATA_(44+I), I =1,4)
     !WRITE(*,*) ( DATA_(44+I), I =1,4)
     ! READ(10,*)  TSTBOT,CSTBOT,TINBOT,CINBOT 
     ! !WRITE(*,*)  TSTBOT,CSTBOT,TINBOT,CINBOT

     !WRITE(*,*)' ISTEADY,ISCHEME,ISLIP,ITBOT,ITTOP,ICBOT,ICTOP '
     READ(10,*)
     
     READ(10,*) ( DATA_(48+I), I =1,7)
     !WRITE(*,*) ( DATA_(48+I), I =1,7)
     ! READ(10,*) ISTEADY,ISCHEME,ISLIP,ITBOT,ITTOP,ICBOT,ICTOP 
     ! !WRITE(*,*) ISTEADY,ISCHEME,ISLIP,ITBOT,ITTOP,ICBOT,ICTOP 

     !WRITE(*,*) ' ITBGRAD,ITTGRAD,ICBGRAD,ICTGRAD '
     READ(10,*)
     READ(10,*) ( DATA_(55+I), I =1,4)
     !WRITE(*,*) ( DATA_(55+I), I =1,4) 
     ! READ(10,*)  ITBGRAD,ITTGRAD,ICBGRAD,ICTGRAD 
     ! !WRITE(*,*)  ITBGRAD,ITTGRAD,ICBGRAD,ICTGRAD 


     !WRITE(*,*)' want to change  x-coordinate (1/0)(y/n) ? '
     READ(10,*)
     
     READ(10,*) ( DATA_(59+I), I =1,1)
     !WRITE(*,*) ( DATA_(59+I), I =1,1)
     ! READ(10,*)  IX
     ! !WRITE(*,*)  IX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IF( IX .EQ. 1 ) THEN  !!!!
!         DX = ARXY/(NX - 1 )!!!!
!         X(1) = 0.          !!!!
!         DO 633 I = 2,NX    !!!!
!            X(I) = X(I-1) + DX!!
! 633     ENDDO                !!
!      ENDIF                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !WRITE(*,*)' want to change  Y-coordinate (1/0)(y/n) ? '
     READ(10,*)
     READ(10,*) ( DATA_(60+I), I =1,1)
     !WRITE(*,*) ( DATA_(60+I), I =1,1)
     ! READ(10,*)  IY
     ! !WRITE(*,*)  IY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !WRITE(*,*)' want to change  Z-coordinate (1/0)(y/n) ? '
     READ(10,*)
     READ(10,*) ( DATA_(61+I), I =1,1)
     !WRITE(*,*) ( DATA_(61+I), I =1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!      IF( IY .EQ. 1 ) THEN
!         DY = 1./(NY - 1 )
!         Y(1) = 0.
!         DO 533 J = 2,NY
!            Y(J) = Y(J-1) + DY
! 533     ENDDO
!      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !WRITE(*,*)' want to see parameters at the terminal( 1/0)(y/n)? '
     READ(10,*)
     READ(10,*) ( DATA_(62+I), I =1,1)
     !WRITE(*,*) ( DATA_(62+I), I =1,1)
     ! READ(10,*)  IPARA
     ! !WRITE(*,*)  IPARA 

     !WRITE(*,*)' want to read data_ from output file(1/0)(y/n)? '
     READ(10,*)
     READ(10,*) ( DATA_(63+I), I =1,1)
     !WRITE(*,*) ( DATA_(63+I), I =1,1)
     ! READ(10,*)  IREAD
     ! !WRITE(*,*)  IREAD

     !WRITE(*,*) ' want to check data_ from output file(1/0)(y/n)? '
     READ(10,*)
     READ(10,*) ( DATA_(64+I), I =1,1)
     !WRITE(*,*) ( DATA_(64+I), I =1,1)
     ! READ(10,*)   ICHECK
     ! !WRITE(*,*)   ICHECK

     !WRITE(*,*)'NUMBER OF TRACER PARTICLES'
     READ(10,*)
     READ(10,*) ( DATA_(65+I), I =1,1)
     !WRITE(*,*) ( DATA_(65+I), I =1,1)
     ! READ(10,*) NTR	
     ! !WRITE(*,*) NTR

     ntr = int(data_(66));
     !WRITE(*,*)' X - COORDINATE FOR TRACER PARTICLES'
     READ(10,*)
     READ(10,*) ( DATA_(66+I), I =1,NTR)
     !WRITE(*,*) ( DATA_(66+I), I =1,NTR)
     ! READ(10,*)(TRX(J),J=1,NTR)
     ! !WRITE(*,*)(TRX(J),J=1,NTR)

     !WRITE(*,*)' Y - COORDINATE FOR TRACER PARTICLES'
     READ(10,*)
     READ(10,*) ( DATA_(66+NTR+I), I =1,NTR)
     !WRITE(*,*) ( DATA_(66+NTR+I), I =1,NTR)
     ! READ(10,*) (TRY(J),J=1,NTR)
     ! !WRITE(*,*) (TRY(J),J=1,NTR)

     !WRITE(*,*) 'NAVTAU,NFINE,nniter,nnfine,NPWD,NXBIG,NYBIG'
     READ(10,*)
     READ(10,*) ( DATA_(66+(2*NTR)+I), I =1,7)
     !WRITE(*,*) ( DATA_(66+(2*NTR)+I), I =1,7)
     ! READ(10,*) NAVTAU,NFINE,NNITER,NNFINE,NPWD,NXBIG,NYBIG
     ! !WRITE(*,*) NAVTAU,NFINE,NNITER,NNFINE,NPWD,NXBIG,NYBIG

     CLOSE(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!READING FINISHES !!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


!     CALL GEOMET(X,Y,NX,NY,XTC,YTC)

     PHIR1 = ( 3.1415926/180)*DATA_(5)
     PHIR2 = ( 3.1415926/180)*DATA_(6)
     DATA_(66 + 2*NTR + 7 + 1) = PHIR1;
     DATA_(66 + 2*NTR + 7 + 2) = PHIR2;
  endif
  
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  
  call MPI_BCAST(DATA_,100,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  ! Get a new communicator for a decomposition of the domain.  Let MPI
  ! find a "good" decomposition
  !
  nx = int(data_(1));ny = int(data_(2));nz = int(data_(3)) ! typecast it into integers
  dims(1) = 0;dims(2) = 0;dims(3) = 0;
  call mpi_dims_create(numprocs, 3, dims, ierr)
  call mpi_cart_create(mpi_comm_world, 3, dims, periods, .true., comm3d, ierr)
  !
  ! Get my position in this communicator
  !
  call mpi_comm_rank(comm3d, myid ,ierr)
  !
  ! My neighbors are now +/- 1 with my rank.  Handle the case of the
  ! boundaries by using MPI_PROCNULL.
  
  !call fnd2dnbrs( comm2d, nbrleft, nbrright, nbrtop, nbrbottom )
  ! Compute the decomposition
  !
!  call fnd2ddecomp( comm2d, nx, ny, sx, ex, sy, ey )
  
  !call mpi_2d(data_,sx,ex,sy,ey,nx,ny,dims,zero_off,comm2d)
  
  call MPI_Comm_free( comm3d, ierr )
  
  t2 = MPI_WTIME();
!  print *, 'total time = ', t2-t1;
  call MPI_FINALIZE(rc)
end program main
