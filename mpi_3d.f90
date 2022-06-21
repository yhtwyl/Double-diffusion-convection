subroutine mpi_3d(data__,sx,ex,sy,ey,sz,ez,nx,ny,nz,dims,zero_off,comm3d)
  use mpi
  implicit none
  integer m,n,p,mp, ivatem
  parameter(m = 102, n = 102, p = 102)
  integer sx, ex, sy, ey, sz, ez, nx, ny, nz, nxyz_, nxyz, n_xyz
  REAL, DIMENSION(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1):: AE,AW,AN,AS,AT,AB,AP
  REAL, DIMENSION(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1):: AVO,AVS,AET,AWT,ANT,AST
  REAL, DIMENSION(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1):: ATOP,ABTM,APT,APOT,APOC
  REAL, DIMENSION(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1):: BE,BW,BN,BS,BT,BB,BP,BVO,BVS,BMAS
  REAL, DIMENSION(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1):: CE,CW,CN,CS,CB,CT,CP,CVO,CVS
  REAL, DIMENSION(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1):: DUM1,DUM2,DUM3,PE,PW,PN,PS,PBTM,PT,PB,pc
  REAL, DIMENSION(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1):: PP,UCP,VCP,WCP,PRS
  REAL, DIMENSION(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1,2):: C,T,UC,VC,WC
  REAL,DIMENSION(M) :: X,XDU,XDV,XDW,XDT,XVEU,XVWU,XWEU,XWWU
  REAL,DIMENSION(M) :: XDEU,XDWU,XDEV,XDWV,XDEW,XDWW,XDET,XDWT,XTC
  REAL,DIMENSION(N) :: Y,YDU,YDV,YDW,YDT,YUNV,YUSV,YDNU,YDSU,YDNW
  REAL,DIMENSION(N) :: YDSW,YDNV,YDSV,YDNT,YDST,YTC,YWNV,YWSV
  REAL,DIMENSION(P) :: Z,ZDU,ZDV,ZDW,ZDT,ZDTT,ZDBT,ZTC,ZDTU
  REAL,DIMENSION(P) :: ZDBU,ZDTV,ZDBV,ZDTW,ZDBW,ZUTW,ZUBW,ZVTW,ZVBW
  REAL,DIMENSION(P) :: TRX,TRY,TRZ
  REAL, DIMENSION(100) :: DATA__
  ! IN FUTURE FOR COMMON DATA
  
  CHARACTER(len=20) temp_fname, conc_fname, uvel_fname, vvel_fname, wvel_fname
!!  CHARACTER(20) FNAME2, FNAME3,FNAME4,FNAME5
  INTEGER ISCHEME,iter, ifine, iflag,itr
  INTEGER LAST_ITER,iSIZE_REC
  integer myid, numprocs, it, rc, comm3d, ierr, stride_EW, stride_NS, stride_TB,k
  integer nbreast, nbrwest, nbrnorth, nbrsouth,nbrtop,nbrbottom,i,ii,ip,itau,j,nitmax,ntau
  integer dims(3), nfine, nnfine, nniter,fh,filetype,gsizes(3),lsizes(3),start_indices(3),memsizes(3),memtype,status; !(MPI_STATUS_SIZE);
  logical periods(3)
  real phir1,umbal,umbal_,vmbal,vmbal_,wmbal,wmbal_,rmmass_,rmmass,enbal_,enbal,conbal_,conbal,tau
  real taumax,td,dtau,dudt_,dudt,dvdt_,dvdt,dwdt_,dwdt
!  external diff2d
  data periods/3*.false./
  integer zero_off

  
  nfine = int(data__(88)); nniter = int(data__(89)); nnfine = int(data__(90));dtau = data__(17);
  taumax = data__(18);iflag = int(data__(19)); ntau = int(data__(26));nitmax = int(data__(27));
  ischeme = int(data__(50))

  ! Get my position in this communicator
  !
  call MPI_COMM_RANK( comm3d, myid, ierr )
  !
  ! My neighbors are now +/- 1 with my rank.  Handle the case of the
  ! boundaries by using MPI_PROCNULL.
  
  call fnd3dnbrs( comm3d, nbreast, nbrwest, nbrnorth, nbrsouth, nbrtop, nbrbottom )
  ! Compute the decomposition
  !
  !call fnd3ddecomp( comm2d, nx, ny, sx, ex, sy, ey )
  ! Create a new, "strided" datatype for the exchange in the "non-contiguous"
  ! direction
  call mpi_Type_vector( (ey-sy+3)*(ez-sz+3), 1, ex-sx+3, mpi_real, stride_EW, ierr )
  call mpi_Type_commit( stride_EW, ierr)
  
  call mpi_Type_vector((ez-sz+3), (ex-sx+3), (ey-sy+3)*(ex-sx+3), mpi_real, stride_NS, ierr)
  call mpi_Type_commit( stride_NS, ierr)

  call mpi_Type_vector((ex-sx+3)*(ey-sy+3), 1 , 1, mpi_real, stride_TB, ierr)
  call mpi_Type_commit( stride_TB, ierr)
  mp = max(m,n,p);
  call Initialise(m,n,p,x,xdu,xdv,xdw,xdt,xveu,xvwu,xweu,xwwu,xdeu,xdwu,xdev,xdwv,xdew,xdww,xdet,xdwt, &
     xtc,y,ydu,ydv,ydw,ydt,yunv,yusv,ywnv,ywsv,ydnu,ydsu,ydnv,ydsv,ydnw,ydsw,ydnt,ydst, &
     ytc,z,zdu,zdv,zdw,zdt,zutw,zubw,zvtw,zvbw,zdtu,zdbu,zdtv,zdbv,zdtw,zdbw,zdtt,zdbt,ztc, &
     nx,ny,nz,sx,ex,sy,ey,sz,ez,data__(20),data__(21),data__(15),data__(16),data__(47),data__(48),&
     t(:,:,:,1),c(:,:,:,1),uc(:,:,:,1),vc(:,:,:,1),wc(:,:,:,1))
  
  ! do k = sz-1,ez+1
  !    do j = sy-1,ey+1
  !       do i = sx-1,ex+1
  !          print *,'before', i,j,k,t(i,j,k,1),myid    
  !       enddo
  !    enddo
  ! enddo

  itau = 0
  tau = 0.
  iter=0
  ifine = 0

90 tau = tau + dtau;
  
! !  print *, 'dtau = ',dtau
  
  itau = itau + 1
  if (itau .eq. nniter) then
     nfine = nnfine;
  endif
  
  if( iflag .eq. 1 ) then
     td =  ((taumax/dtau) + 0.000001) 
     ii = ifix(td)      
     if( itau .gt. ii ) then
        tau = taumax
        itau = itau - 1
!        print *, 'dtau',dtau
        go to 250
     endif
  endif

  itr = 0;
140 iter = iter + 1;
  itr = itr + 1;

  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  
  
  call exchng3(uc(:,:,:,1), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  call exchng3(vc(:,:,:,1), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth &
       , nbrbottom,nbrtop)
  call exchng3(wc(:,:,:,1), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  call exchng3(T(:,:,:,1), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth, &
       nbrbottom,nbrtop)
  call exchng3(C(:,:,:,1), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth, &
       nbrbottom,nbrtop)
  
  call ucoef1(m,n,p,uc(:,:,:,1),uc(:,:,:,2),vc(:,:,:,2),wc(:,:,:,2),t(:,:,:,2),c(:,:,:,2),nx,ny,nz,sx,ex,sy,ey,sz,ez,data__(13),&
       data__(11),data__(9),data__(17),data__(39),data__(50),data__(6),data__(7),data__(8),xdu,xdeu,xdwu,xveu,xvwu,xweu,xwwu,ydu,&
       ydnu,ydsu,zdu,zdtu,zdbu,ae,aw,an,as,at,ab,ap,avo,avs,myid)

  call vcoef1(m,n,p,uc(:,:,:,2),vc(:,:,:,1),vc(:,:,:,2),wc(:,:,:,2),t(:,:,:,2),c(:,:,:,2),nx,ny,nz,sx,ex,sy,ey,sz,ez,data__(13),&
     data__(11),data__(9),data__(17),data__(39),data__(50),data__(6),data__(7),data__(8),xdv,xdev,xdwv,ydv,ydnv, &
     ydsv,yusv,yunv,ywsv,ywnv,zdu,zdtv,zdbv,be,bw,bn,bs,bt,bb,bp,bvo,bvs,myid)

  call wcoef1(m,n,p,uc(:,:,:,2),vc(:,:,:,2),wc(:,:,:,1),wc(:,:,:,2),t(:,:,:,2),c(:,:,:,2),nx,ny,nz,sx,ex,sy,ey,sz,ez,data__(13),&
     data__(11),data__(9),data__(17),data__(39),data__(50),data__(6),data__(7),data__(8),xdw,xdew,xdww,ydw,ydnw, &
     ydsw,zdw,zdtw,zdbw,zubw,zutw,zvbw,zvtw,ce,cw,cn,cs,ct,cb,cp,cvo,cvs,myid)
  
  ! do k = sz-1,ez+1
  !    do j = sy-1,ey+1
  !       do i = sx-1,ex+1
  !          print *,'Before', i,j,k,bp(i,j,k),myid
  !       enddo
  !    enddo
  ! enddo
  
  call exchng3(ap, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  
  call exchng3(bp, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  
  call exchng3(cp, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)


  call ucap(m,n,p,uc(:,:,:,2),nx,ny,nz,sx,ex,sy,ey,sz,ez,ae,aw,an,as,at,ab,ap,avo,avs,ucp,myid)

  call vcap(m,n,p,vc(:,:,:,2),nx,ny,nz,sx,ex,sy,ey,sz,ez,be,bw,bn,bs,bt,bb,bp,bvo,bvs,vcp,myid)

  call wcap(m,n,p,wc(:,:,:,2),nx,ny,nz,sx,ex,sy,ey,sz,ez,ce,cw,cn,cs,ct,cb,cp,cvo,cvs,wcp,myid)

  call exchng3(ucp, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)

  call exchng3(vcp, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  
  call exchng3(wcp, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  
  call prcoff(m,n,p,ap,bp,cp,ucp,vcp,wcp,pe,pw,pn,ps,pt,pbtm,pp,pb,xdt,ydt,zdt,nx,ny,nz,sx,ex,sy,ey,sz,ez,myid)
  
  do ii = 1,3
     call presrx(mp,prs,pe,pw,pn,ps,pt,pbtm,pp,pb,nx,ny,nz,sx,ex,sy,ey,sz,ez)
     call exchng3(prs, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
          nbrbottom,nbrtop)
     
     call presry(mp,prs,pe,pw,pn,ps,pt,pbtm,pp,pb,nx,ny,nz,sx,ex,sy,ey,sz,ez)
     call exchng3(prs, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
          nbrbottom,nbrtop)
     
     call presrz(mp,prs,pe,pw,pn,ps,pt,pbtm,pp,pb,nx,ny,nz,sx,ex,sy,ey,sz,ez)
     call exchng3(prs, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  enddo
  
  ! do k = sz,ez
  !    do j = sy,ey
  !       do i = sx,ex
  !          print *,'Before', i,j,k,pe(i,j,k),pw(i,j,k),pn(i,j,k),ps(i,j,k),pt(i,j,k),pbtm(i,j,k),myid
  !       enddo
  !    enddo
  ! enddo
  
  call umomex(m,n,p,mp,uc(:,:,:,2),prs,aw,ae,as,an,at,ab,ap,avo,avs,ydu,zdu,nx,ny,nz,sx,ex,sy,ey,sz,ez)
  call umomey(m,n,p,mp,uc(:,:,:,2),prs,aw,ae,as,an,at,ab,ap,avo,avs,ydu,zdu,nx,ny,nz,sx,ex,sy,ey,sz,ez)
  call umomez(m,n,p,mp,uc(:,:,:,2),prs,aw,ae,as,an,at,ab,ap,avo,avs,ydu,zdu,nx,ny,nz,sx,ex,sy,ey,sz,ez,umbal_,comm3d,myid)

  call exchng3(uc(:,:,:,2), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  call mpi_allreduce(umbal_, umbal, 1, mpi_real, mpi_sum, comm3d, ierr)
  if (myid .eq. 0) umbal = sqrt(umbal/((nz-1)*(ny-1)*(nx-2))) ! ITS TOO HIGH
  call mpi_bcast(umbal,1,mpi_real,0,comm3d,ierr)
!  print *, '192_mpi_3d = ',umbal_,umbal

  call vmomex(m,n,p,mp,vc(:,:,:,2),prs,bw,be,bs,bn,bt,bb,bp,bvo,bvs,xdv,zdv,nx,ny,nz,sx,ex,sy,ey,sz,ez)
  call vmomey(m,n,p,mp,vc(:,:,:,2),prs,bw,be,bs,bn,bt,bb,bp,bvo,bvs,xdv,zdv,nx,ny,nz,sx,ex,sy,ey,sz,ez)
  call vmomez(m,n,p,mp,vc(:,:,:,2),prs,bw,be,bs,bn,bt,bb,bp,bvo,bvs,xdv,zdv,nx,ny,nz,sx,ex,sy,ey,sz,ez,vmbal_,comm3d,myid)

  call exchng3(vc(:,:,:,2), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  call mpi_allreduce(vmbal_, vmbal, 1, mpi_real, mpi_sum, comm3d, ierr)
  if (myid .eq. 0) vmbal = sqrt(vmbal/((nz-1)*(ny-2)*(nx-1))) ! ITS TOO HIGH
  call mpi_bcast(vmbal,1,mpi_real,0,comm3d,ierr)
  !  print *, '211_mpi_3d = ',vmbal_,vmbal
  
  call wmomex(m,n,p,mp,wc(:,:,:,2),prs,cw,ce,cs,cn,ct,cb,cp,cvo,cvs,xdw,ydw,nx,ny,nz,sx,ex,sy,ey,sz,ez)
  call wmomey(m,n,p,mp,wc(:,:,:,2),prs,cw,ce,cs,cn,ct,cb,cp,cvo,cvs,xdw,ydw,nx,ny,nz,sx,ex,sy,ey,sz,ez)
  call wmomez(m,n,p,mp,wc(:,:,:,2),prs,cw,ce,cs,cn,ct,cb,cp,cvo,cvs,xdw,ydw,nx,ny,nz,sx,ex,sy,ey,sz,ez,wmbal_,comm3d,myid)

  call exchng3(wc(:,:,:,2), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  call mpi_allreduce(wmbal_, wmbal, 1, mpi_real, mpi_sum, comm3d, ierr)
  if (myid .eq. 0) wmbal = sqrt(wmbal/((nz-2)*(ny-1)*(nx-1))) ! ITS TOO HIGH
  call mpi_bcast(wmbal,1,mpi_real,0,comm3d,ierr)
!  print *, '222_mpi_3d = ',wmbal_,wmbal
!   call mpi_allreduce(vmbal_, vmbal, 1, mpi_real, mpi_sum, comm2d, ierr)
!   if (myid .eq. 0) vmbal = sqrt(vmbal/((ny-2)*(nx-1))) ! ITS TOO HIGH
!   call mpi_bcast(vmbal,1,mpi_real,0,comm2d,ierr)

  call mass(m,n,p,uc(:,:,:,2),vc(:,:,:,2),wc(:,:,:,2),nx,ny,nz,sx,ex,sy,ey,sz,ez,bmas,iter,data__(25),data__(26),rmmass_,&
       xdt,ydt,zdt,nxyz_)
  call mpi_allreduce(nxyz_, n_xyz, 1, mpi_int, mpi_sum,comm3d, ierr)
  call mpi_allreduce(rmmass_, rmmass, 1, mpi_real, mpi_sum, comm3d, ierr)

  if (myid .eq. 0) then
     nxyz = (nx - 1)*(ny -1)*(nz -1) - n_xyz;
     if (nxyz .ne. 0) rmmass = sqrt(rmmass/nxyz)
!     print *, rmmass
  endif
  call mpi_bcast(rmmass,1,mpi_real,0,comm3d,ierr)

  call prcorx(mp,pe,pw,pn,ps,pt,pbtm,pc,pp,sx,ex,sy,ey,sz,ez,nx,ny,nz,bmas) 
  call prcory(mp,pe,pw,pn,ps,pt,pbtm,pc,pp,sx,ex,sy,ey,sz,ez,nx,ny,nz,bmas)
  call prcorz(mp,pe,pw,pn,ps,pt,pbtm,pc,pp,sx,ex,sy,ey,sz,ez,nx,ny,nz,bmas)
  call exchng3(pc, sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  
  ! do k = sz,ez
  !    do j = sy,ey
  !       do i = sx,ex
  !          print *,'Before', i,j,k,pc(i,j,k),myid
  !       enddo
  !    enddo
  ! enddo
  
  call uvel(m,n,p,uc(:,:,:,2),pc,nx,ny,nz,sx,ex,sy,ey,sz,ez,ap,ydu,zdu)
  call vvel(m,n,p,vc(:,:,:,2),pc,nx,ny,nz,sx,ex,sy,ey,sz,ez,bp,xdv,zdv)
  call wvel(m,n,p,wc(:,:,:,2),pc,nx,ny,nz,sx,ex,sy,ey,sz,ez,cp,xdw,ydw)
  call exchng3(uc(:,:,:,2), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  call exchng3(vc(:,:,:,2), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  call exchng3(wc(:,:,:,2), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)

  call species_coeff(m,n,p,uc(:,:,:,2),vc(:,:,:,2),wc(:,:,:,2),t(:,:,:,1),aet,awt,ant,ast,atop,abtm,apt,apot,xdet,&
       xdt,xdwt,ydt,ydnt,ydst,zdt,zdbt,zdtt,nx,ny,nz,sx,ex,sy,ey,sz,ez,data__(13),data__(14),data__(17),data__(50),&
       data__(56),data__(57),'t_')
  call species_calc_x(mp,t(:,:,:,2),aet,awt,ant,ast,atop,abtm,apt,apot,sx,ex,sy,ey,sz,ez,nx,ny,nz,data__(43),data__(19),&
       data__(52),data__(53),data__(56),data__(57))!! check for ytc and ztc
  call species_calc_y(mp,t(:,:,:,2),aet,awt,ant,ast,atop,abtm,apt,apot,sx,ex,sy,ey,sz,ez,nx,ny,nz,data__(43),data__(19),&
       data__(52),data__(53),data__(56),data__(57))!! check for ytc and ztc
  call species_calc_z(mp,t(:,:,:,2),aet,awt,ant,ast,atop,abtm,apt,apot,sx,ex,sy,ey,sz,ez,nx,ny,nz,data__(43),data__(19),&
       data__(52),data__(53),data__(56),data__(57),enbal_)!! check for ytc and ztc
  call exchng3(t(:,:,:,2), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)

  call mpi_allreduce(enbal_, enbal, 1, mpi_real, mpi_sum, comm3d, ierr)
 
  if (myid .eq. 0)  enbal = sqrt(enbal/((nx-1)*(ny-1)*(nz-1)))
!  if (myid .eq. 0) print *, 'enbal', enbal
  call mpi_bcast(enbal,1,mpi_real,0,comm3d,ierr)

  call species_coeff(m,n,p,uc(:,:,:,2),vc(:,:,:,2),wc(:,:,:,2),c(:,:,:,1),aet,awt,ant,ast,atop,abtm,apt,apot,xdet,xdt,&
       xdwt,ydt,ydnt,ydst,zdt,zdbt,zdtt,nx,ny,nz,sx,ex,sy,ey,sz,ez,data__(13),data__(14),data__(17),data__(50),data__(58),&
       data__(59),'c_')
  call species_calc_x(mp,c(:,:,:,2),aet,awt,ant,ast,atop,abtm,apt,apot,sx,ex,sy,ey,sz,ez,nx,ny,nz,data__(44),data__(19),&
       data__(54),data__(55),data__(58),data__(59))!! check for ytc and ztc
  call species_calc_y(mp,c(:,:,:,2),aet,awt,ant,ast,atop,abtm,apt,apot,sx,ex,sy,ey,sz,ez,nx,ny,nz,data__(44),data__(19),&
       data__(54),data__(55),data__(58),data__(59))!! check for ytc and ztc
  call species_calc_z(mp,c(:,:,:,2),aet,awt,ant,ast,atop,abtm,apt,apot,sx,ex,sy,ey,sz,ez,nx,ny,nz,data__(44),data__(19),&
       data__(54),data__(55),data__(58),data__(59),conbal_)!! check for ytc and ztc
  call exchng3(t(:,:,:,2), sx, ex, sy, ey, sz, ez, comm3d, stride_EW, stride_NS, stride_TB, nbrwest, nbreast, nbrsouth, nbrnorth,&
       nbrbottom,nbrtop)
  call mpi_allreduce(conbal_, conbal, 1, mpi_real, mpi_sum, comm3d, ierr)
  
  if (myid .eq. 0)  conbal = sqrt(conbal/((nx-1)*(ny-1)*(nz-1)))
  !  if (myid .eq. 0) print *, 'conbal',conbal
  call mpi_bcast(conbal,1,mpi_real,0,comm3d,ierr)
  
  do k = sz,ez
     do j=sy,ey
        do i=sx,ex
           pc(i,j,k)=0.0
        enddo
     enddo
  enddo

  if (myid .eq. 0) then
     if(int(data__(63)) .eq. 1) then
        call screen(data__(9),data__(11),tau,data__(17),iter,itr,itau,rmmass,nxyz,enbal,conbal,umbal,vmbal,wmbal)
     endif
  endif

  if( iflag .eq. 0 ) then
    ! if( iter .gt. itmass ) then  
    !    if( rmmass .le. ermass .or. nxy .le. nacell ) then
    !       go to 210
    !    endif
    ! endif
  else
     if( itau .le. ntau ) then
        if( itr .eq. nitmax ) go to 210
     else
        if( (itr .eq. nitmax) .or. (rmmass .le. data__(33)) )go to 210
     endif
  endif

  go to 140
  
210 continue
 
  call ududt(m,n,p,uc(:,:,:,2),prs,nx,ny,nz,sx,ex,sy,ey,sz,ez,ae,aw,an,as,at,ab,avs,xdu,ydu,zdu,dudt_)
  call mpi_allreduce(dudt_, dudt, 1, mpi_real, mpi_sum, comm3d, ierr)
  if (myid .eq. 0) dudt = sqrt(dudt/((nx-2)*(ny-1)*(nz-1)))
  call mpi_bcast(dudt,1,mpi_real,0,comm3d,ierr)
  call vdvdt(m,n,p,vc(:,:,:,2),prs,nx,ny,nz,sx,ex,sy,ey,sz,ez,be,bw,bn,bs,bt,bb,bvs,xdv,ydv,zdv,dvdt_)
  call mpi_allreduce(dvdt_, dvdt, 1, mpi_real, mpi_sum, comm3d, ierr)
  if (myid .eq. 0) dvdt = sqrt(dvdt/((nx-1)*(ny-2)*(nz-1)))
  call mpi_bcast(dvdt,1,mpi_real,0,comm3d,ierr)
  call wdwdt(m,n,p,wc(:,:,:,2),prs,nx,ny,nz,sx,ex,sy,ey,sz,ez,ce,cw,cn,cs,ct,cb,cvs,xdw,ydw,zdw,dwdt_)
  call mpi_allreduce(dwdt_, dwdt, 1, mpi_real, mpi_sum, comm3d, ierr)
  if (myid .eq. 0) dwdt = sqrt(dwdt/((nx-1)*(ny-1)*(nz-2)))
  call mpi_bcast(dwdt,1,mpi_real,0,comm3d,ierr)
  
!   ! print *,'sx,ex,sy,ey', sx ,ex, sy, ey
!  print *, 'dudt,dvdt',dudt,dvdt,dwdt
  do k = sz,ez
     do j = sy,ey
        do i = sx,ex
           uc(i,j,k,1) = uc(i,j,k,2)
           vc(i,j,k,1) = vc(i,j,k,2)
           wc(i,j,k,1) = wc(i,j,k,2)
           t(i,j,k,1) = t(i,j,k,2)
           c(i,j,k,1) = c(i,j,k,2)
!           print *, i,j,k, t(i,j,k,1),t(i,j,k,2)
        enddo
     enddo
  enddo
print *,'nfine = ', nfine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!  FILE WRITING !!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (itau .eq. (nfine*(itau/nfine))) then
     call filename('t', itau, temp_fname)
     call filename('c', itau, conc_fname)
     call filename('u', itau, uvel_fname)
     call filename('v', itau, vvel_fname)
     call filename('w', itau, wvel_fname)
     gsizes(1) = nx;gsizes(2) = ny;gsizes(3) = nz; ! no of rows and columns in global array
     lsizes(1) = ex - sx + 1; lsizes(2) = ey - sy + 1; lsizes(3) = ez - sz + 1;
     start_indices(1) = sx-1;start_indices(2) = sy-1;start_indices(3) = sz-1;
     call mpi_type_create_subarray(3,gsizes,lsizes,start_indices, mpi_order_fortran, mpi_real, filetype, ierr);
     call mpi_type_commit(filetype, ierr);
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call mpi_file_open(mpi_comm_world, temp_fname,mpi_mode_create + mpi_mode_wronly &
          , mpi_info_null, fh, ierr);! mpi_mode_create .or. /home/fluids1/Documents/MPI/DDC/
     call mpi_file_set_view(fh, 0_mpi_offset_kind , mpi_real, filetype, "native", mpi_info_null, ierr); ! 0_mpi_offset_kind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !    call write2file(nx,ny,sx,sy,ex,ey,t(:,:,2),fh,my_id)
     !    call write2file(nx,ny,sx,sy,ex,ey,t(:,:,2),fh,)
     memsizes(1) = ex-sx+1; memsizes(2) = ey-sy+1; memsizes(3) = ez-sz+1;
     start_indices(1) = 0; start_indices(2) = 0; start_indices(3) = 0;
     !    print *,"mym = ",myid,memsizes(1),memsizes(2),lsizes(1),lsizes(2),start_indices(1),start_indices(2)
     call mpi_type_create_subarray(3, memsizes, lsizes, start_indices, mpi_order_fortran, mpi_real, memtype, ierr);
     call mpi_type_commit(memtype,ierr);
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call mpi_file_write_all(fh, t(sx:ex,sy:ey,sz:ez,1), 1, memtype, status, ierr)
     call mpi_file_close(fh,ierr);
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call mpi_file_open(mpi_comm_world, conc_fname,mpi_mode_create + mpi_mode_wronly &
          , mpi_info_null, fh, ierr);! mpi_mode_create .or. /home/fluids1/Documents/MPI/DDC/
     call mpi_file_set_view(fh, 0_mpi_offset_kind , mpi_real, filetype, "native", mpi_info_null, ierr); ! 0_mpi_offset_kind


     call mpi_file_write_all(fh, c(sx:ex,sy:ey,sz:ez,1), 1, memtype, status, ierr)
     call mpi_file_close(fh,ierr);
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call mpi_file_open(mpi_comm_world, uvel_fname,mpi_mode_create + mpi_mode_wronly &
          , mpi_info_null, fh, ierr);! mpi_mode_create .or. /home/fluids1/Documents/MPI/DDC/
     call mpi_file_set_view(fh, 0_mpi_offset_kind , mpi_real, filetype, "native", mpi_info_null, ierr); ! 0_mpi_offset_kind


     call mpi_file_write_all(fh, uc(sx:ex,sy:ey,sz:ez,1), 1, memtype, status, ierr)
     call mpi_file_close(fh,ierr);
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call mpi_file_open(mpi_comm_world, vvel_fname,mpi_mode_create + mpi_mode_wronly &
          , mpi_info_null, fh, ierr);! mpi_mode_create .or. /home/fluids1/Documents/MPI/DDC/
     call mpi_file_set_view(fh, 0_mpi_offset_kind , mpi_real, filetype, "native", mpi_info_null, ierr); ! 0_mpi_offset_kind


     call mpi_file_write_all(fh, vc(sx:ex,sy:ey,sz:ez,1), 1, memtype, status, ierr)
     call mpi_file_close(fh,ierr);
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call mpi_file_open(mpi_comm_world, wvel_fname,mpi_mode_create + mpi_mode_wronly &
          , mpi_info_null, fh, ierr);! mpi_mode_create .or. /home/fluids1/Documents/MPI/DDC/
     call mpi_file_set_view(fh, 0_mpi_offset_kind , mpi_real, filetype, "native", mpi_info_null, ierr); ! 0_mpi_offset_kind


     call mpi_file_write_all(fh, wc(sx:ex,sy:ey,sz:ez,1), 1, memtype, status, ierr)
     call mpi_file_close(fh,ierr);
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     call mpi_type_free(filetype,ierr)
     call mpi_type_free(memtype,ierr)
 endif
  
  go to 90
  
250 ip = 0 

  if((dudt .le. data__(36)) .and. (dvdt .le. data__(37)) .and. (dwdt .le. data__(38)))then
     if (myid .eq. 0) then
        write(*,*)'solution is converged for transient eqn. after'
        write(*,*)'time step = ',tau
     endif
  else
     if (myid .eq. 0) then
        write(*,*)'solution is not converged for transient eqn. after'
        write(*,*)'time step = ',tau
     endif
  endif

260 continue

  call mpi_type_free(stride_EW, ierr)
  call mpi_type_free(stride_NS, ierr)
  call mpi_type_free(stride_TB, ierr)
  return
end subroutine mpi_3d




