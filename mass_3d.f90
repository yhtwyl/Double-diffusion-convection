subroutine mass(m,n,p,uc_,vc_,wc_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,bmas_,iter_,itmass_,conmas_,rmmass__,xdt_,ydt_,zdt_,nxyz__)
!  use mpi
!  parameter(m=22,n=22)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,nxyz__,iter_,ierr_,myid_,itmass_
  real rmmass__,aa_,a1_,a2_,a3_,a4_,a5_,a6_,conmas_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_,vc_,wc_,bmas_,bmsin_
  real, dimension(m) :: xdt_
  real, dimension(n) :: ydt_
  real, dimension(p) :: zdt_
!  nxy = (nx_-1)*(ny_-1);
  nxyz__ = 0
  aa_ = 0.0
  do k = sz_,ez_
     if ((k .gt. 1) .and. (k .lt. nz_+1)) then
        do j = sy_,ey_
           if ((j .gt. 1) .and. (j .lt. ny_+1)) then
              do i = sx_,ex_
                 if ((i .gt. 1) .and. (i .lt. nx_+1)) then
                    bmas_(i,j,k) = (uc_(i-1,j,k) - uc_(i,j,k))*ydt_(j)*zdt_(k) + (vc_(i,j-1,k)-vc_(i,j,k))*xdt_(i)*zdt_(k) +&
                         (wc_(i,j,k-1)-wc_(i,j,k))*xdt_(i)*ydt_(j);
                    a1_ = uc_(i-1,j,k)*ydt_(j)*zdt_(k);a2_ = uc_(i,j,k)*ydt_(j)*zdt_(k);
                    a3_ = vc_(i,j-1,k)*xdt_(i)*zdt_(k);a4_ = vc_(i,j,k)*xdt_(i)*zdt_(k);
                    a5_ = wc_(i,j,k-1)*xdt_(i)*ydt_(j);a6_ = wc_(i,j,k)*xdt_(i)*ydt_(j);
                    bmsin_(i,j,k)=MAX(aa_,a1_) + MAX(aa_,-a2_) + MAX(aa_,a3_) + MAX(aa_,-a4_) + max(aa_,a5_) + max(aa_,-a6_);
                    if (abs(bmsin_(i,j,k)) .lt. conmas_) nxyz__ = nxyz__ + 1
                    !              print *,'mass' ,i,j,a1_,a2_,a3_,a4_
                 endif
              enddo
           endif
        enddo
     endif
  enddo

  if (iter_ .gt. itmass_) then
     rmmass__ = 0.0
     do k = sz_,ez_
        if ((k .gt. 1) .and. (k .lt. nz_+1)) then
           do j = sy_,ey_
              if ((j .gt. 1) .and. (j .lt. ny_+1)) then
                 do i = sx_,ex_
                    if ((i .gt. 1) .and. (i .lt. nx_+1)) then
                       if (abs(bmsin_(i,j,k)) .gt. conmas_) then
                          rmmass__ = rmmass__ + (bmas_(i,j,k)/bmsin_(i,j,k))**2;
                       endif
                    endif
                 enddo
              endif
           enddo
        endif
     enddo
  endif
end subroutine mass
