subroutine uvel(m,n,p,uc_,pc_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,ap_,ydu_,zdu_)
!  parameter(m=22,n=22)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_,pc_,ap_
!  real, dimension(m) :: a__,b__,c__,d__,x__
  real, dimension(n) :: ydu_
  real, dimension(p) :: zdu_
  
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           if (j .gt. 1) then
              do i = sx_,ex_
                 if ((i .gt. 1) .and. (i .lt. nx_))then
                    uc_(i,j,k) = uc_(i,j,k) + (ydu_(j)*zdu_(k)/ap_(i,j,k))*(pc_(i,j,k)- pc_(i+1,j,k))
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine uvel

subroutine vvel(m,n,p,vc_,pc_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,bp_,xdv_,zdv_)
!  parameter(m=22,n=22)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: vc_,pc_,bp_
  real, dimension(m) :: xdv_
  real, dimension(p) :: zdv_
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           if ((j .gt. 1) .and. (j .lt. ny_)) then
              do i = sx_,ex_
                 if (i .gt. 1) then
                    vc_(i,j,k) = vc_(i,j,k) + (xdv_(i)*zdv_(k)/bp_(i,j,k))*(pc_(i,j,k) - pc_(i,j+1,k));
!              print * ,i,j,vc_(i,j)
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine vvel

subroutine wvel(m,n,p,wc_,pc_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,cp_,xdw_,ydw_)
!  parameter(m=22,n=22)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: wc_,pc_,cp_
  real, dimension(m) :: xdw_
  real, dimension(n) :: ydw_
  do k = sz_,ez_
     if ((k .gt. 1) .and. (k .lt. nz_)) then
        do j = sy_,ey_
           if (j .gt. 1) then
              do i = sx_,ex_
                 if (i .gt. 1) then
                    wc_(i,j,k) = wc_(i,j,k) + (xdw_(i)*ydw_(j)/cp_(i,j,k))*(pc_(i,j,k) - pc_(i,j,k+1));
!              print * ,i,j,vc_(i,j)
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine wvel
