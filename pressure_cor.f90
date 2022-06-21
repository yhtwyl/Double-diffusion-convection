subroutine prcorx(mp,pe_,pw_,pn_,ps_,pt_,pbtm_,pc_,pp_,sx_,ex_,sy_,ey_,sz_,ez_,nx_,ny_,nz_,bmas_) 
!  parameter(m=22,n=22)
  integer mp,i,j,k,i_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  real, dimensIon(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: pe_,pw_,pn_,ps_,pt_,pbtm_,pc_,pp_,bmas_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  start_ = sx_ - 2;
  do k = sz_,ez_
     if (k .gt. 1) then
        do j=sy_,ey_
           if (j .gt. 1) then
              i_ = 0
              do i=sx_-1,ex_+1
                 if ((i .gt. 1) .and. (i .lt. nx_+1)) then
                    i_ = i_ + 1;
                    if ((i .ne. sx_-1) .and. (i .ne. ex_+1)) then
                       d__(i_) = pn_(i,j,k)*pc_(i,j+1,k) + ps_(i,j,k)*pc_(i,j-1,k) + pt_(i,j,k)*pc_(i,j,k+1) + &
                            pbtm_(i,j,k)*pc_(i,j,k-1) + bmas_(i,j,k);
                       b__(i_) = pp_(i,j,k);
                       if (i .eq. 2) then
                          a__(i_) = 0.0;c__(i_) = -pe_(i,j,k);
                       elseif(i .eq. nx_) then
                          a__(i_) = -pw_(i,j,k);c__(i_) = 0.0;
                       else
                          a__(i_) = -pw_(i,j,k);c__(i_) = -pe_(i,j,k);
                       endif
                    else
                       a__(i_) = 0.0;b__(i_) = 1.0;c__(i_)=0.0;d__(i_)=pc_(i,j,k);
                    endif
                    if (i .eq. 2) start_ = 1;
                 endif
              enddo
              call tridag(a__,b__,c__,d__,i_,x__,mp,'prcorx')
              do i = 1,i_
                 pc_(start_ + i,j,k) = x__(i)
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine prcorx

subroutine prcory(mp,pe_,pw_,pn_,ps_,pt_,pbtm_,pc_,pp_,sx_,ex_,sy_,ey_,sz_,ez_,nx_,ny_,nz_,bmas_)
!  parameter(m=22,n=22)
  integer mp,i,j,k,j_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: pe_,pw_,pn_,ps_,pt_,pbtm_,pc_,pp_,bmas_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  start_ = sy_ - 2;
  do k = sz_,ez_
     if (k .gt. 1) then
        do i=sx_,ex_
           if (i .gt. 1) then
              j_ = 0
              do j=sy_-1,ey_+1
                 if ((j .gt. 1) .and. (j .lt. ny_+1)) then
                    j_ = j_ + 1;
                    if ((j .ne. sy_-1) .and. (j .ne. ey_+1)) then
                       d__(j_) = pe_(i,j,k)*pc_(i+1,j,k) + pw_(i,j,k)*pc_(i-1,j,k) + pt_(i,j,k)*pc_(i,j,k+1) + &
                            pbtm_(i,j,k)*pc_(i,j,k-1) + bmas_(i,j,k);
                       b__(j_) = pp_(i,j,k);
                       if (j .eq. 2) then
                          a__(j_) = 0.0;c__(j_) = -pn_(i,j,k);
                       elseif(j .eq. ny_) then
                          a__(j_) = -ps_(i,j,k);c__(j_) = 0.0;
                       else
                          a__(j_) = -ps_(i,j,k);c__(j_) = -pn_(i,j,k);
                       endif
                    else
                       a__(j_) = 0.0;b__(j_) = 1.0;c__(j_)=0.0;d__(j_)=pc_(i,j,k);
                    endif
                    if (j .eq. 2) start_ = 1;
                 endif
              enddo
              call tridag(a__,b__,c__,d__,j_,x__,mp,'prcory')
              do j = 1,j_
                 pc_(i,start_ + j,k) = x__(j)
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine prcory

subroutine prcorz(mp,pe_,pw_,pn_,ps_,pt_,pbtm_,pc_,pp_,sx_,ex_,sy_,ey_,sz_,ez_,nx_,ny_,nz_,bmas_)
!  parameter(m=22,n=22)
  integer mp,i,j,k,j_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,start_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: pe_,pw_,pn_,ps_,pt_,pbtm_,pc_,pp_,bmas_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  start_ = sz_ - 2;
  do j=sy_,ey_
     if (j .gt. 1) then
        do i=sx_,ex_
           if (i .gt. 1) then
              k_ = 0
              do k = sz_-1,ez_+1
                 if ((k .gt. 1) .and. (k .lt. nz_+1)) then
                    k_ = k_ + 1;
                    if ((k .ne. sz_-1) .and. (k .ne. ez_+1)) then
                       d__(k_) = pe_(i,j,k)*pc_(i+1,j,k) + pw_(i,j,k)*pc_(i-1,j,k) + pn_(i,j,k)*pc_(i,j+1,k) + &
                            ps_(i,j,k)*pc_(i,j-1,k) + bmas_(i,j,k);
                       b__(k_) = pp_(i,j,k);
                       if (k .eq. 2) then
                          a__(k_) = 0.0;c__(k_) = -pt_(i,j,k);
                       elseif(k .eq. nz_) then
                          a__(k_) = -pbtm_(i,j,k);c__(k_) = 0.0;
                       else
                          a__(k_) = -pbtm_(i,j,k);c__(k_) = -pt_(i,j,k);
                       endif
                    else
                       a__(k_) = 0.0;b__(k_) = 1.0;c__(k_)=0.0;d__(k_)=pc_(i,j,k);
                    endif
                    if (k .eq. 2) start_ = 1;
                 endif
              enddo
              call tridag(a__,b__,c__,d__,k_,x__,mp,'prcorz')
              do k = 1,k_
                 pc_(i,j,start_ + k) = x__(k)
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine prcorz
