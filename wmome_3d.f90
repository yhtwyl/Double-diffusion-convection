subroutine wmomex(m,n,p,mp,wc_,prs_,cw_,ce_,cs_,cn_,ct_,cb_,cp_,cvo_,cvs_,xdw_,ydw_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_)
  integer m,n,p,mp,i,j,k,i_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: wc_,prs_,ce_,cw_,cn_,cs_,cb_,ct_,cp_,cvo_,cvs_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  real, dimension(m) :: xdw_
  real, dimension(n) :: ydw_
  start_ = sx_-2
  do k = sz_,ez_
     if ((k .gt. 1) .and. (k .lt. nz_)) then
        do j = sy_,ey_
           if (j .gt. 1) then
              i_ = 0
              do i = sx_-1,ex_+1
                 if ((i .gt. 1) .and. (i .lt. nx_+1)) then
                    i_ = i_ + 1;
                    if ((i .ne. sx_-1) .and. (i .ne. ex_+1)) then
                       d__(i_) = cn_(i,j,k)*wc_(i,j+1,k) + cs_(i,j,k)*wc_(i,j-1,k) + ct_(i,j,k)*wc_(i,j,k+1) + cb_(i,j,k) &
                            *wc_(i,j,k-1)+ ( prs_(i,j,k) - prs_(i,j,k+1))*xdw_(i)*ydw_(j) + cvo_(i,j,k) + cvs_(i,j,k);
                       b__(i_)=  cp_(i,j,k)
                       a__(i_)= -cw_(i,j,k)
                       c__(i_)= -ce_(i,j,k)
                    elseif ((i .eq. sx_-1) .or. (i .eq. ex_+1)) then
                       a__(i_) = 0.0;b__(i_) = 1.0;c__(i_)=0.0;d__(i_)=wc_(i,j,k);
                    endif
                 endif
                 if (i .eq. 2) start_ = 1
              enddo
              call tridag(a__,b__,c__,d__,i_,x__,mp,'wmomex')
              do i=1,i_
                 wc_(start_ + i, j, k) = x__(i);
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine wmomex

subroutine wmomey(m,n,p,mp,wc_,prs_,cw_,ce_,cs_,cn_,ct_,cb_,cp_,cvo_,cvs_,xdw_,ydw_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_)
  !  Parameter(m=22,n=22)
  integer m,n,p,mp,i,j,k,j_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: wc_,prs_,ce_,cw_,cn_,cs_,ct_,cb_,cp_,cvo_,cvs_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  real, dimension(m) :: xdw_
  real, dimension(n) :: ydw_
  start_ = sy_ - 2;
  do k = sz_,ez_
     if ((k .gt. 1) .and. (k .lt. nz_)) then
        do i = sx_,ex_
           if (i .gt. 1) then
              j_ = 0;
              do j = sy_-1,ey_+1
                 if ((j .gt. 1) .and. (j .lt. ny_ + 1)) then
                    j_ = j_ + 1;
                    if ((j .ne. sy_-1) .and. (j .ne. ey_+1)) then
                       d__(j_) = ce_(i,j,k)*wc_(i+1,j,k) + cw_(i,j,k)*wc_(i-1,j,k) + ct_(i,j,k)*wc_(i,j,k+1) + cb_(i,j,k) &
                            *wc_(i,j,k-1)+ (prs_(i,j,k) - prs_(i,j,k+1))*xdw_(i)*ydw_(j) + cvo_(i,j,k) + cvs_(i,j,k);
                       b__(j_) = cp_(i,j,k);a__(j_) = -cs_(i,j,k);c__(j_)=-cn_(i,j,k);
                    elseif ((j .eq. sy_-1) .or. (j .eq. ey_+1)) then
                       a__(j_) = 0.0;b__(j_) = 1.0;c__(j_)=0.0;d__(j_)=wc_(i,j,k);
                    endif
                 endif
                 if (j .eq. 2) start_ = 1;
                 !              print *, 'umomex',i,j,a__(i_),b__(i_),c__(i_),d__(i_)
              enddo
              call tridag(a__,b__,c__,d__,j_,x__,mp,'wmomey')
              do j=1,j_
                 wc_(i,start_ + j,k) = x__(j);
              !           print *,'creator prs_',i,start_+i,j,prs_(start_+i,j)
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine wmomey


subroutine wmomez(m,n,p,mp,wc_,prs_,cw_,ce_,cs_,cn_,ct_,cb_,cp_,cvo_,cvs_,xdw_,ydw_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,&
     wmbal__,comm3d,myid_)
  integer m,n,p,mp,i,j,k,k_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_,ierr_,myid_,comm3d
  real wmbal_,wmbal__,wbal_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: wc_,prs_,ce_,cw_,cn_,cs_,cb_,ct_,cp_,cvo_,cvs_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  real, dimension(m) :: xdw_
  real, dimension(n) :: ydw_
  start_ = sz_-2
  do j = sy_,ey_
     if (j .gt. 1) then
        do i = sx_,ex_
           if (i .gt. 1) then
              k_ = 0
              do k = sz_-1,ez_+1
                 if ((k .gt. 0) .and. (k .lt. nz_+1)) then
                    k_ = k_ + 1;
                    if ((k .gt. 1) .and. (k .lt. nz_)) then
                       if ((k .ne. sz_ - 1) .and. (k .ne. ez_ + 1)) then
                          d__(k_) = ce_(i,j,k)*wc_(i+1,j,k) + cw_(i,j,k)*wc_(i-1,j,k) + cn_(i,j,k)*wc_(i,j+1,k) + &
                               cs_(i,j,k)*wc_(i,j-1,k) + ( prs_(i,j,k) - prs_(i,j,k+1))*xdw_(i)*ydw_(j) + cvo_(i,j,k) + cvs_(i,j,k);
                          b__(k_)=  cp_(i,j,k)
                          a__(k_)= -cb_(i,j,k)
                          c__(k_)= -ct_(i,j,k)
                       elseif ((k .eq. (sz_-1)) .or. (k .eq. (ez_+1))) then
                          a__(k_) = 0.0;b__(k_) = 1.0;c__(k_)=0.0;d__(k_)=wc_(i,j,k);
                       endif
                    elseif ((k .eq. 1) .or. (k .eq. nz_)) then
                       a__(k_) = 0.0;b__(k_) = 1.0;c__(k_)=0.0;d__(k_)=0.0;
                    endif
                    if (k .eq. 1) start_ = 0;
                 endif
              enddo
              call tridag(a__,b__,c__,d__,k_,x__,mp,'vmomez')
              do k = 1,k_
                 wc_(i,j,start_ + k) = x__(k);
              enddo
           endif
        enddo
     endif
  enddo
  !! to find x_momentum residuals
  wmbal__ = 0.0
  do k = sz_,ez_
     if ((k .gt. 1) .and. (k .lt. nz_)) then
        do j =  sy_,ey_
           if (j .gt. 1) then
              do i =  sx_,ex_
                 if (i .gt. 1) then
                    wbal_ = cp_(i,j,k)*wc_(i,j,k) - ce_(i,j,k)*wc_(i+1,j,k) - cw_(i,j,k)*wc_(i-1,j,k) - cn_(i,j,k)*wc_(i,j+1,k) -  &
                         cs_(i,j,k)*wc_(i,j-1,k) - ct_(i,j,k)*wc_(i,j,k+1) - cb_(i,j,k)*wc_(i,j,k-1) - cvo_(i,j,k) - cvs_(i,j,k)&
                         - (prs_(i,j,k) - prs_(i,j,k+1))*xdw_(i)*ydw_(j);
                    wmbal__ = wmbal__ + (wbal_**2);
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine wmomez
