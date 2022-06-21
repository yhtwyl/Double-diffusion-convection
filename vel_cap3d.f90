subroutine ucap(m,n,p,uc_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,ae_,aw_,an_,as_,at_,ab_,ap_,avo_,avs_,ucp_,my_id)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,my_id
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: ae_,aw_,an_,as_,at_,ab_,ap_,avo_,avs_,uc_,ucp_

  do k = sz_,ez_
     if (k .gt. 1) then
        do j= sy_,ey_
           if (j .gt. 1) then
              do i = sx_,ex_
                 if ((i .gt. 1) .and. (i .lt. nx_)) then
                    ucp_(i,j,k)=(ae_(i,j,k)*uc_(i+1,j,k) + aw_(i,j,k)*uc_(i-1,j,k) + an_(i,j,k)*uc_(i,j+1,k) + as_(i,j,k)*&
                         uc_(i,j-1,k) + at_(i,j,k)*uc_(i,j,k+1) + ab_(i,j,k)*uc_(i,j,k-1)+avo_(i,j,k) + avs_(i,j,k))/&
                         ap_(i,j,k)
                 endif
              enddo
           endif
        enddo
     endif
  enddo
end subroutine ucap


subroutine vcap(m,n,p,vc_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,be_,bw_,bn_,bs_,bt_,bb_,bp_,bvo_,bvs_,vcp_,my_id)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,my_id
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: be_,bw_,bn_,bs_,bt_,bb_,bp_,bvo_,bvs_,vc_,vcp_

  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           if ((j .gt. 1) .and. (j .lt. ny_)) then
              do  i = sx_,ex_
                 if (i .gt. 1) then
                    vcp_(i,j,k) = (be_(i,j,k)*vc_(i+1,j,k) + bw_(i,j,k)*vc_(i-1,j,k) + bn_(i,j,k)*vc_(i,j+1,k) + &
                         bs_(i,j,k)*vc_(i,j-1,k) + bt_(i,j,k)*vc_(i,j,k+1) + bb_(i,j,k)*vc_(i,j,k-1) + bvo_(i,j,k)&
                         + bvs_(i,j,k))/bp_(i,j,k)
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine vcap

subroutine wcap(m,n,p,wc_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,ce_,cw_,cn_,cs_,ct_,cb_,cp_,cvo_,cvs_,wcp_,my_id)
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,my_id
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: ce_,cw_,cn_,cs_,ct_,cb_,cp_,cvo_,cvs_,wc_,wcp_

  do k = sz_,ez_
     if ((k .gt. 1) .and. (k .lt. nz_)) then
        do j = sy_,ey_
           if (j .gt. 1) then
              do  i = sx_,ex_
                 if (i .gt. 1) then
                    wcp_(i,j,k) = (ce_(i,j,k)*wc_(i+1,j,k) + cw_(i,j,k)*wc_(i-1,j,k) + cn_(i,j,k)*wc_(i,j+1,k) + &
                         cs_(i,j,k)*wc_(i,j-1,k) + ct_(i,j,k)*wc_(i,j,k+1) + cb_(i,j,k)*wc_(i,j,k-1) + cvo_(i,j,k)&
                         + cvs_(i,j,k))/cp_(i,j,k)
!                    print *, i,j,k,cvo_(i,j,k),cvs_(i,j,k),my_id
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine wcap
