subroutine vmomex(m,n,p,mp,vc_,prs_,bw_,be_,bs_,bn_,bt_,bb_,bp_,bvo_,bvs_,xdv_,zdv_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_)
  integer m,n,p,mp,i,j,k,i_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: vc_,prs_,be_,bw_,bn_,bs_,bb_,bt_,bp_,bvo_,bvs_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  real, dimension(m) :: xdv_
  real, dimension(p) :: zdv_
  start_ = sx_-2
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           if ((j .gt. 1) .and. (j .lt. ny_)) then
              i_ = 0
              do i = sx_-1,ex_+1
                 if ((i .gt. 1) .and. (i .lt. nx_+1)) then
                    i_ = i_ + 1;
                    if ((i .ne. sx_-1) .and. (i .ne. ex_+1)) then
                       d__(i_) = bn_(i,j,k)*vc_(i,j+1,k) + bs_(i,j,k)*vc_(i,j-1,k) + bt_(i,j,k)*vc_(i,j,k+1) + bb_(i,j,k) &
                            *vc_(i,j,k-1)+ ( prs_(i,j,k) - prs_(i,j+1,k))*xdv_(i)*zdv_(k) + bvo_(i,j,k) + bvs_(i,j,k);
                       b__(i_)=  bp_(i,j,k)
                       a__(i_)= -bw_(i,j,k)
                       c__(i_)= -be_(i,j,k)
                    elseif ((i .eq. sx_-1) .or. (i .eq. ex_+1)) then
                       a__(i_) = 0.0;b__(i_) = 1.0;c__(i_)=0.0;d__(i_)=vc_(i,j,k);
                    endif
                 endif
                 if (i .eq. 2) start_ = 1
              enddo
              call tridag(a__,b__,c__,d__,i_,x__,mp,'vmomex')
              do i=1,i_
                 vc_(start_ + i, j, k) = x__(i);
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine vmomex

subroutine vmomey(m,n,p,mp,vc_,prs_,bw_,be_,bs_,bn_,bt_,bb_,bp_,bvo_,bvs_,xdv_,zdv_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_)
!  Parameter(m=22,n=22)
  integer m,n,p,mp,i,j,k,j_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: vc_,prs_,be_,bw_,bn_,bs_,bt_,bb_,bp_,bvo_,bvs_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  real, dimension(m) :: xdv_
  real, dimension(p) :: zdv_
  start_ = sy_ - 2;
  do k = sz_,ez_
     if (k .gt. 1) then
        do i = sx_,ex_
           if (i .gt. 1) then
              j_ = 0;
              do j = sy_-1,ey_+1
                 if ((j .gt. 0) .and. (j .lt. ny_ + 1)) then
                    j_ = j_ + 1;
                    if ((j .gt. 1) .and. (j .lt. ny_)) then
                       if ((j .ne. sy_-1) .and. (j .ne. ey_+1)) then
                          d__(j_) = be_(i,j,k)*vc_(i+1,j,k) + bw_(i,j,k)*vc_(i-1,j,k) + bt_(i,j,k)*vc_(i,j,k+1) + bb_(i,j,k) &
                               *vc_(i,j,k-1)+ (prs_(i,j,k) - prs_(i,j+1,k))*xdv_(i)*zdv_(k) + bvo_(i,j,k) + bvs_(i,j,k);
                          b__(j_) = bp_(i,j,k);a__(j_) = -bs_(i,j,k);c__(j_)=-bn_(i,j,k);
                       elseif ((j .eq. sy_-1) .or. (j .eq. ey_+1)) then
                          a__(j_) = 0.0;b__(j_) = 1.0;c__(j_)=0.0;d__(j_)=vc_(i,j,k);
                       endif
                    elseif ((j .eq. 1) .or. (j .eq. ny_)) then
                       a__(j_) = 0.0;b__(j_) = 1.0;c__(j_)=0.0;d__(j_)=0.0;
                    endif
                    if (j .eq. 1) start_ = 0;
                    !              print *, 'umomex',i,j,a__(i_),b__(i_),c__(i_),d__(i_)
                 endif
              enddo
              call tridag(a__,b__,c__,d__,j_,x__,mp,'vmomey')
              do j=1,j_
                 vc_(i,start_ + j,k) = x__(j);
                 !           print *,'creator prs_',i,start_+i,j,prs_(start_+i,j)
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine vmomey


subroutine vmomez(m,n,p,mp,vc_,prs_,bw_,be_,bs_,bn_,bt_,bb_,bp_,bvo_,bvs_,xdv_,zdv_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,vmbal__,&
     comm3d,myid_)
  integer m,n,p,mp,i,j,k,k_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_,ierr_,myid_,comm3d
  real vmbal_,vmbal__,vbal_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: vc_,prs_,be_,bw_,bn_,bs_,bb_,bt_,bp_,bvo_,bvs_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  real, dimension(m) :: xdv_
  real, dimension(p) :: zdv_
  start_ = sz_-2
  do j = sy_,ey_
     if ((j .gt. 1) .and. (j .lt. ny_)) then
        do i = sx_,ex_
           if (i .gt. 1) then
              k_ = 0
              do k = sz_-1,ez_+1
                 if ((k .gt. 1) .and. (k .lt. nz_+1)) then
                    k_ = k_ + 1;
                    if ((k .ne. sz_-1) .and. (k .ne. ez_+1)) then
                       d__(k_) = be_(i,j,k)*vc_(i+1,j,k) + bw_(i,j,k)*vc_(i-1,j,k) + bn_(i,j,k)*vc_(i,j+1,k) + &
                            bs_(i,j,k)*vc_(i,j-1,k) + ( prs_(i,j,k) - prs_(i,j+1,k))*xdv_(i)*zdv_(k) + bvo_(i,j,k)&
                            + bvs_(i,j,k);
                       b__(k_)=  bp_(i,j,k)
                       a__(k_)= -bb_(i,j,k)
                       c__(k_)= -bt_(i,j,k)
                    elseif ((k .eq. (sz_-1)) .or. (k .eq. (ez_+1))) then
                       a__(k_) = 0.0;b__(k_) = 1.0;c__(k_)=0.0;d__(k_)=vc_(i,j,k);
                    endif
                 endif
                 if (k .eq. 2) start_ = 1
              enddo
              call tridag(a__,b__,c__,d__,k_,x__,mp,'vmomez')
              do k = 1,k_
                 vc_(i,j,start_ + k) = x__(k);
              enddo
           endif
        enddo
     endif
  enddo
  !! to find x_momentum residuals
  vmbal__ = 0.0
  do k = sz_,ez_
     if (k .gt. 1) then
        do j =  sy_,ey_
           if ((j .gt. 1) .and. (j .lt. ny_)) then
              do i =  sx_,ex_
                 if (i .gt. 1) then
                    vbal_ = bp_(i,j,k)*vc_(i,j,k) - be_(i,j,k)*vc_(i+1,j,k) - bw_(i,j,k)*vc_(i-1,j,k) - bn_(i,j,k)*vc_(i,j+1,k) -  &
                         bs_(i,j,k)*vc_(i,j-1,k) - bt_(i,j,k)*vc_(i,j,k+1) - bb_(i,j,k)*vc_(i,j,k-1) - bvo_(i,j,k) - bvs_(i,j,k)&
                         - (prs_(i,j,k) - prs_(i,j+1,k))*xdv_(i)*zdv_(k);
                    vmbal__ = vmbal__ + (vbal_**2);
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine vmomez
