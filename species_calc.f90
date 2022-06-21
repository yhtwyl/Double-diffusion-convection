! C*****************************************************************************
! C   SUBROUTINE TO CALCULATE TEMPERATURE FIELD
! C--> ( this routine calculates temp. and concentration field for all 
! C-->   types of bouadary conditions ( fixed/non-zero grad./zero grad.)
! C-->   at bottom and top )
! C-->  IBGRAD = 0  bottom at zero gradient condition
! C-->         = 1  bottom at fixed/non-zero gradient condition( same as
! C-->                                              ITBGRAD/ICBGRAD )
! C-->  ITGRAD = 0  top at zero gradient condition
! C-->         = 1  top at fixed/non-zerogradient condition( same as  
! C-->                                              ITTGRAD/ICTGRAD )    
! C-->  IBOT = 1  bottom at fixed boundary condition
! C-->       = 0  bottom at gradient boundary condition( zero or non-zero)
! C-->  ITOP = 1  top at fixed boundary condition
! C-->       = 0  top at gradient boundary condition( zero or non-zero )
! C*******************************************************************************

subroutine species_calc_x(mp,t_,aet_,awt_,ant_,ast_,atop_,abtm_,apt_,apot_,sx_,ex_,sy_,ey_,sz_,ez_,nx_,ny_,nz_,relax_,iflag_,&
     ibot_,itop_,ibgrad_,itgrad_)
!  parameter(m=22,n=22)
  integer m,n,p,mp,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,i_,start_
  integer ibot_,itop_,ibgrad_,itgrad_,iflag_
  real relax_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: t_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: aet_,awt_,ant_,ast_,atop_,abtm_,apt_,apot_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  start_ = sx_ - 2
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           i_ = 0
           if (j .gt. 1) then
              do i = sx_-1,ex_+1
                 if ((i .gt. 1) .and. (i .lt. nx_+1))then
                    i_ = i_ + 1;
                    if ((i .ne. (sx_ - 1)) .and. (i .ne. (ex_ + 1))) then
                       a__(i_) = -awt_(i,j,k)
                       b__(i_) =  apt_(i,j,k) 
                       c__(i_) = -aet_(i,j,k)
                       d__(i_) =  ant_(i,j,k)*t_(i,j+1,k) + ast_(i,j,k)*t_(i,j-1,k) + atop_(i,j,k)*t_(i,j,k+1) + &
                            abtm_(i,j,k)*t_(i,j,k-1) + apot_(i,j,k)
                    else
                       a__(i_) = 0.0;b__(i_) = 1.0;c__(i_)=0.0;d__(i_)=t_(i,j,k);
                    endif
                    if(i .eq. 2) start_ = 1;
                    !              print *, i,j, aet_(i,j),awt_(i,j),ant_(i,j),ast_(i,j),apot_(i,j)
                 endif
              enddo
              call tridag(a__,b__,c__,d__,i_,x__,mp,'species_calc_x')
              do i=1,i_
                 t_(start_ + i,j,k) = x__(i);
              enddo
           endif
        enddo
     endif
  enddo
!!!!!!!!!!!! no flux at boundary !!!!!!!!
!!!! check it for speed of code
  if (sx_ .eq. 1 ) then
     t_(1,sy_:ey_,sz_:ez_) = t_(2,sy_:ey_,sz_:ez_);
  elseif (ex_ .eq. nx_) then
     t_(nx_+1,sy_:ey_,sz_:ez_) = t_(nx_,sy_:ey_,sz_:ez_);
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  return
end subroutine species_calc_x



subroutine species_calc_y(mp,t_,aet_,awt_,ant_,ast_,atop_,abtm_,apt_,apot_,sx_,ex_,sy_,ey_,sz_,ez_,nx_,ny_,nz_,relax_,iflag_,ibot_,&
     itop_,ibgrad_,itgrad_)
!  parameter(m=22,n=22)
  integer mp,i,j,k,j_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  integer ibot_,itop_,ibgrad_,itgrad_,iflag_
  real ebal_,enbal_,enbal__,relax_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: t_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: aet_,awt_,ant_,ast_,atop_,abtm_,apt_,apot_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  start_ = sy_ -2
  do k = sz_,ez_
     if (k .gt. 1) then
        do i = sx_,ex_
           j_ = 0
           if (i .gt. 1) then
              do j = sy_-1,ey_+1
                 if ((j .gt. 1) .and. (j .lt. ny_+1))then
                    j_ = j_ + 1;
                    if ((j .ne. (sy_ - 1)) .and. (j .ne. (ey_ + 1))) then
                       a__(j_) = -ast_(i,j,k)
                       b__(j_) =  apt_(i,j,k) 
                       c__(j_) = -ant_(i,j,k)
                       d__(j_) =  aet_(i,j,k)*t_(i+1,j,k) + awt_(i,j,k)*t_(i-1,j,k) + atop_(i,j,k)*t_(i,j,k+1) + &
                            abtm_(i,j,k)*t_(i,j,k-1) + apot_(i,j,k)
                    else
                       a__(j_) = 0.0;b__(j_) = 1.0;c__(j_)=0.0;d__(j_)=t_(i,j,k);
                    endif
                    if (j .eq. 2) start_ = 1;
                 endif
              enddo
              call tridag(a__,b__,c__,d__,j_,x__,mp,'species_calc_y')
              do j=1,j_
                 t_(i,start_ + j,k) = x__(j);
           !           print *,'creator prs_',i,start_+i,j,prs_(start_+i,j)
              enddo
           endif
        enddo
     endif
  enddo
!!!!!!!!!!!! no flux at boundary !!!!!!!!
!!!! check it for speed of code
  if (sy_ .eq. 1) then
     t_(sx_:ex_,1,sz_:ez_) = t_(sx_:ex_,2,sz_:ez_);
  elseif (sy_ .eq. ny_) then
     t_(sx_:ex_,ny_+1,sz_:ez_) = t_(sx_:ex_,ny_,sz_:ez_);
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  return
end subroutine species_calc_y

subroutine species_calc_z(mp,t_,aet_,awt_,ant_,ast_,atop_,abtm_,apt_,apot_,sx_,ex_,sy_,ey_,sz_,ez_,nx_,ny_,nz_,relax_,iflag_,ibot_,&
     itop_,ibgrad_,itgrad_,enbal__)
  integer mp,i,j,k,k_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  integer ibot_,itop_,ibgrad_,itgrad_,iflag_
  real ebal_,enbal__,relax_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: t_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: aet_,awt_,ant_,ast_,atop_,abtm_,apt_,apot_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  start_ = sz_ -2
  do j = sy_,ey_
     if (j .gt. 1) then
        do i = sx_,ex_
           if (i .gt. 1) then
              k_ = 0;
              do k = sz_-1,ez_+1
                 if ((k .gt. 1) .and. (k .lt. nz_+1))then
                    k_ = k_ + 1;
                    if ((k .ne. (sz_ - 1)) .and. (k .ne. (ez_ + 1))) then
                       a__(k_) = -ast_(i,j,k)
                       b__(k_) =  apt_(i,j,k) 
                       c__(k_) = -ant_(i,j,k)
                       d__(k_) =  aet_(i,j,k)*t_(i+1,j,k) + awt_(i,j,k)*t_(i-1,j,k) + atop_(i,j,k)*t_(i,j,k+1) + &
                            abtm_(i,j,k)*t_(i,j,k-1) + apot_(i,j,k)
                    else
                       a__(k_) = 0.0;b__(k_) = 1.0;c__(k_)=0.0;d__(k_)=t_(i,j,k);
                    endif
                    if (k .eq. 2) start_ = 1;
                 endif
              enddo
              call tridag(a__,b__,c__,d__,k_,x__,mp,'species_calc_z')
              do k=1,k_
                 t_(i,j,start_ + k) = x__(k);
!                 print *,'creator prs_',i,start_+i,j,prs_(start_+i,j)
              enddo
           endif
        enddo
     endif
  enddo
!!!!!!!!!!!! no flux at boundary !!!!!!!!
  !!!! check it for speed of code
  if (sz_ .eq. 1) then
     t_(sx_:ex_,sy_:ey_,1) = t_(sx_:ex_,sy_:ey_,2);
  elseif (sz_ .eq. nz_) then
     t_(sx_:ex_,sy_:ey_,nz_+1) = t_(sx_:ex_,sy_:ey_,nz_);
  endif
  !!!!!!!!!!!! no flux at boundary !!!!!!!!
  
!!! to find energy residuals !!!!!!!!!!!!!!!
  enbal__ = 0.0;
  do k = sz_,ez_
     do j = sy_,ey_
        do i = sx_,ex_
           ebal_ = apt_(i,j,k)*t_(i,j,k) - aet_(i,j,k)*t_(i+1,j,k) - awt_(i,j,k)*t_(i-1,j,k) - ant_(i,j,k)*t_(i,j+1,k)&
                - ast_(i,j,k)*t_(i,j-1,k) - atop_(i,j,k)*t_(i,j,k+1) - abtm_(i,j,k)*t_(i,j,k-1) - apot_(i,j,k);
           enbal__ = enbal__ + (ebal_**2);
        enddo
     enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  return
end subroutine species_calc_z
