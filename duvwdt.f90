subroutine ududt(m,n,p,uc_,prs_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,ae_,aw_,an_,as_,at_,ab_,avs_,xdu_,ydu_,zdu_,dudt__)
!  parameter(m = 22, n = 22)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_
  real balu_,dudt__
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_,prs_,ae_,aw_,an_,as_,at_,ab_,avs_
  real, dimension(m) :: xdu_
  real, dimension(n) :: ydu_
  real, dimension(p) :: zdu_
  dudt__ = 0.0
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           if (j .gt. 1) then
              do i = sx_,ex_
                 if ((i .gt. 1 ) .and. (i .lt. nx_)) then
                    balu_ = aw_(i,j,k)*(uc_(i-1,j,k)-uc_(i,j,k)) - ae_(i,j,k)*(uc_(i,j,k)-uc_(i+1,j,k)) + as_(i,j,k)*(uc_(i,j-1,k)-&
                         uc_(i,j,k)) - an_(i,j,k)*(uc_(i,j,k)-uc_(i,j+1,k)) + ab_(i,j,k)*(uc_(i,j,k-1) - uc_(i,j,k)) - &
                         at_(i,j,k)*(uc_(i,j,k) - uc_(i,j,k+1)) + (prs_(i,j,k)-prs_(i+1,j,k))*ydu_(j)*zdu_(k) +avs_(i,j,k);
                    balu_ = balu_/(xdu_(i)*ydu_(j)*zdu_(k));
                    dudt__ = dudt__ + (balu_**2);
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine ududt

subroutine vdvdt(m,n,p,vc_,prs_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,be_,bw_,bn_,bs_,bt_,bb_,bvs_,xdv_,ydv_,zdv_,dvdt__)
!  parameter(m = 22, n = 22)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_
  real balv_,dvdt__
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: vc_,prs_,be_,bw_,bn_,bs_,bt_,bb_,bvs_
  real, dimension(m) :: xdv_
  real, dimension(n) :: ydv_
  real, dimension(p) :: zdv_
  dvdt__ = 0.0
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           if ((j .gt. 1) .and. (j .lt. ny_)) then
              do i = sx_,ex_
                 if (i .gt. 1 ) then
                    balv_ = bw_(i,j,k)*(vc_(i-1,j,k)-vc_(i,j,k)) - be_(i,j,k)*(vc_(i,j,k)-vc_(i+1,j,k)) + bs_(i,j,k)*(vc_(i,j-1,k)-&
                         vc_(i,j,k)) - bn_(i,j,k)*(vc_(i,j,k)-vc_(i,j+1,k)) + bb_(i,j,k)*(vc_(i,j,k-1) - vc_(i,j,k)) - &
                         bt_(i,j,k)*(vc_(i,j,k) - vc_(i,j,k+1)) + (prs_(i,j,k)-prs_(i,j+1,k))*xdv_(i)*zdv_(k) + bvs_(i,j,k);
                    balv_ = balv_/(xdv_(i)*ydv_(j)*zdv_(k));
                    dvdt__ = dvdt__ + (balv_**2);
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine vdvdt

subroutine wdwdt(m,n,p,wc_,prs_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,ce_,cw_,cn_,cs_,ct_,cb_,cvs_,xdw_,ydw_,zdw_,dwdt__)
!  parameter(m = 22, n = 22)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_
  real balw_,dwdt__
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: wc_,prs_,ce_,cw_,cn_,cs_,ct_,cb_,cvs_
  real, dimension(m) :: xdw_
  real, dimension(n) :: ydw_
  real, dimension(p) :: zdw_
  dvdt__ = 0.0
  do k = sz_,ez_
     if ((k .gt. 1) .and. (k .lt. nz_)) then
        do j = sy_,ey_
           if (j .gt. 1) then
              do i = sx_,ex_
                 if (i .gt. 1 ) then
                    balw_ = cw_(i,j,k)*(wc_(i-1,j,k)-wc_(i,j,k)) - ce_(i,j,k)*(wc_(i,j,k)-wc_(i+1,j,k)) + cs_(i,j,k)*(wc_(i,j-1,k)-&
                         wc_(i,j,k)) - cn_(i,j,k)*(wc_(i,j,k)-wc_(i,j+1,k)) + cb_(i,j,k)*(wc_(i,j,k-1) - wc_(i,j,k)) - &
                         ct_(i,j,k)*(wc_(i,j,k) - wc_(i,j,k+1)) + (prs_(i,j,k)-prs_(i,j,k+1))*xdw_(i)*ydw_(j) + cvs_(i,j,k);
                    balw_ = balw_/(xdw_(i)*ydw_(j)*zdw_(k));
                    dwdt__ = dwdt__ + (balw_**2);
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine wdwdt
