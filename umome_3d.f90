subroutine umomex(m,n,p,mp,uc_,prs_,aw_,ae_,as_,an_,at_,ab_,ap_,avo_,avs_,ydu_,zdu_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_)
!  Parameter(m=22,n=22)
  integer m,n,p,mp,i,j,k,i_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_,prs_,ae_,aw_,an_,as_,at_,ab_,ap_,avo_,avs_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  real, dimension(n) :: ydu_
  real, dimension(p) :: zdu_
  start_ = sx_ - 2;
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           if (j .gt. 1) then
              i_ = 0
              do i = sx_-1,ex_-1
                 if ((i .gt. 0) .and. (i .lt. nx_ + 1)) then
                    i_ = i_ + 1;
                    if ((i .gt. 1) .and. (i .lt. nx_)) then
                       if ((i .ne. sx_-1) .and. (i .ne. ex_+1)) then
                          d__(i_) = an_(i,j,k)*uc_(i,j+1,k) + as_(i,j,k)*uc_(i,j-1,k) + at_(i,j,k)*uc_(i,j,k+1) + ab_(i,j,k) &
                               *uc_(i,j,k-1)+ (prs_(i,j,k) - prs_(i+1,j,k))*ydu_(j)*zdu_(k) + avo_(i,j,k) + avs_(i,j,k);
                          b__(i_) = ap_(i,j,k);a__(i_) = -aw_(i,j,k);c__(i_)=-ae_(i,j,k);
                       elseif ((i .eq. sx_-1) .or. (i .eq. ex_+1)) then
                          a__(i_) = 0.0;b__(i_) = 1.0;c__(i_)=0.0;d__(i_)=uc_(i,j,k);
                       endif
                    elseif ((i .eq. 1) .or. (i .eq. nx_)) then
                       a__(i_) = 0.0;b__(i_) = 1.0;c__(i_)=0.0;d__(i_)=0.0;
                    endif
                    if (i .eq. 1) start_ = 0;
                    !              print *, 'umomex',i,j,a__(i_),b__(i_),c__(i_),d__(i_)
                 endif
              enddo
              call tridag(a__,b__,c__,d__,i_,x__,mp,'umomex')
              do i=1,i_
                 uc_(start_ + i,j,k) = x__(i);
                 !           print *,'creator prs_',i,start_+i,j,prs_(start_+i,j)
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine umomex 

subroutine umomey(m,n,p,mp,uc_,prs_,aw_,ae_,as_,an_,at_,ab_,ap_,avo_,avs_,ydu_,zdu_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_)
  integer m,n,p,mp,i,j,k,j_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_,prs_,ae_,aw_,an_,as_,ab_,at_,ap_,avo_,avs_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  real, dimension(n) :: ydu_
  real, dimension(p) :: zdu_
  start_ = sy_-2
  do k = sz_,ez_
     if (k .gt. 1) then
        do i = sx_,ex_
           if (i .gt. 1 .and. i .lt. nx_) then
              j_ = 0
              do j = sy_-1,ey_+1
                 if ((j .gt. 1) .and. (j .lt. ny_+1)) then
                    j_ = j_ + 1;
                    if ((j .ne. sy_-1) .and. (j .ne. ey_+1)) then
                       d__(j_) = ae_(i,j,k)*uc_(i+1,j,k) + aw_(i,j,k)*uc_(i-1,j,k) + at_(i,j,k)*uc_(i,j,k+1) + ab_(i,j,k) &
                            *uc_(i,j,k-1)+ ( prs_(i,j,k) - prs_(i+1,j,k))*ydu_(j)*zdu_(k) + &
                            avo_(i,j,k) + avs_(i,j,k);
                       b__(j_)=  ap_(i,j,k)
                       a__(j_)= -as_(i,j,k)
                       c__(j_)= -an_(i,j,k)
                    elseif ((j .eq. sy_-1) .or. (j .eq. ey_+1)) then
                       a__(j_) = 0.0;b__(j_) = 1.0;c__(j_)=0.0;d__(j_)=uc_(i,j,k);
                    endif
                 endif
                 if (j .eq. 2) start_ = 1
              enddo
              call tridag(a__,b__,c__,d__,j_,x__,mp,'umomey')
              do j=1,j_
                 uc_(i,start_ + j,k) = x__(j);
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine umomey

subroutine umomez(m,n,p,mp,uc_,prs_,aw_,ae_,as_,an_,at_,ab_,ap_,avo_,avs_,ydu_,zdu_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,&
     umbal__,comm3d,myid_)
  integer m,n,p,mp,i,j,k,k_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_,ierr_,myid_,comm3d
  real umbal__,ubal_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_,prs_,ae_,aw_,an_,as_,ab_,at_,ap_,avo_,avs_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  real, dimension(n) :: ydu_
  real, dimension(p) :: zdu_
  start_ = sz_-2
  do j = sy_,ey_
     if (j .gt. 1) then
        do i = sx_,ex_
           if ((i .gt. 1) .and. (i .lt. nx_)) then
              k_ = 0
              do k = sz_-1,ez_+1
                 if ((k .gt. 1) .and. (k .lt. nz_+1)) then
                    k_ = k_ + 1;
                    if ((k .ne. sz_-1) .and. (k .ne. ez_+1)) then
                       d__(k_) = ae_(i,j,k)*uc_(i+1,j,k) + aw_(i,j,k)*uc_(i-1,j,k) + an_(i,j,k)*uc_(i,j+1,k) +  &
                            as_(i,j,k)*uc_(i,j-1,k) + ( prs_(i,j,k) - prs_(i+1,j,k))*ydu_(j)*zdu_(k) + &
                            avo_(i,j,k) + avs_(i,j,k);
                       b__(k_)=  ap_(i,j,k)
                       a__(k_)= -ab_(i,j,k)
                       c__(k_)= -at_(i,j,k)
                    elseif ((k .eq. sz_-1) .or. (k .eq. ez_+1)) then
                       a__(k_) = 0.0;b__(k_) = 1.0;c__(k_)=0.0;d__(k_)=uc_(i,j,k);
                    endif
                 endif
                 if (k .eq. 2) start_ = 1
              enddo
              call tridag(a__,b__,c__,d__,k_,x__,mp,'umomez')
              do k=1,k_
                 uc_(i,j,start_ + k) = x__(k);
              enddo
           endif
        enddo
     endif
  enddo
!! to find x_momentum residuals
  umbal__ = 0.0
  do k = sz_,ez_
     if (k .gt. 1) then
        do j =  sy_,ey_
           if (j .gt. 1) then
              do i =  sx_,ex_
                 if (i .gt. 1 .and. i .lt. nx_) then
                    ubal_ = ap_(i,j,k)*uc_(i,j,k) - ae_(i,j,k)*uc_(i+1,j,k) - aw_(i,j,k)*uc_(i-1,j,k) - an_(i,j,k)*uc_(i,j+1,k) -  &
                         as_(i,j,k)*uc_(i,j-1,k) - at_(i,j,k)*uc_(i,j,k+1) - ab_(i,j,k)*uc_(i,j,k-1) - avo_(i,j,k) - avs_(i,j,k) - &
                         (prs_(i,j,k) - prs_(i+1,j,k))*ydu_(j)*zdu_(k);
                    umbal__ = umbal__ + (ubal_**2);
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine umomez

