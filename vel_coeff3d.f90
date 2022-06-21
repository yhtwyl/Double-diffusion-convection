subroutine ucoef1(m,n,p,uc_p,uc_,vc_,wc_,t_,c_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,prndl_,grasm_,grast_,dtau_, &
     sign_,ischeme_,phir1_,phir2_,phir3_,xdu_,xdeu_,xdwu_,xveu_,xvwu_,xweu_,xwwu_,ydu_,ydnu_, &
     ydsu_,zdu_,zdtu_,zdbu_,ae_,aw_,an_,as_,at_,ab_,ap_,avo_,avs_,my_id)

  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,ischeme_,my_id
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_p,uc_,vc_,wc_,t_,c_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: ae_,aw_,an_,as_,at_,ab_,ap_,avo_,avs_
  real prndl_,grasm_,grast_,dtau_,phir1_,phir2_,phir3_,de_,ue_,fe_,pe_,aa_,aut_,au_,gamma_
  real dw_,uw_,fw_,pw_,dn_,vn_,fn_,pn_,ds_,vs_,fs_,ps_,dt_,wt_,ft_,pt_,db_,wb_,fb_,pb_,st_,sc_,st1_,sc1_,sign_
  real, dimension(m) :: xdu_,xdeu_,xdwu_,xveu_,xvwu_,xweu_,xwwu_
  real, dimension(n) :: ydu_,ydnu_,ydsu_
  real, dimension(p) :: zdu_,zdtu_,zdbu_
  gamma_= prndl_

  do k = sz_,ez_
     if (k .gt. 1) then
        do j=sy_,ey_
           if (j .gt. 1) then
              do i=sx_,ex_
                 if ((i .gt. 1) .and. (i .lt. nx_)) then
                    de_=(gamma_*ydu_(j)*zdu_(k))/xdeu_(i)
                    ue_=(uc_(i,j,k) + uc_(i+1,j,k))*0.5
                    fe_=ue_*ydu_(j)*zdu_(k)
                    pe_=fe_/de_
                    aa_=0.
                    if( ischeme_ .eq. 1 ) then 
                       aut_= (1. - 0.1*abs(pe_))**5
                    else
                       aut_= 1 - 0.5*abs(pe_)
                    endif
                    au_=amax1(aa_,aut_)
                    ae_(i,j,k)=de_*au_+amax1(-fe_,aa_)
                    !c-->
                    dw_=(gamma_*ydu_(j)*zdu_(k))/xdwu_(i)
                    uw_=(uc_(i,j,k)+uc_(i-1,j,k))*0.5
                    fw_=uw_*ydu_(j)*zdu_(k)
                    pw_=fw_/dw_
                    if( ischeme_ .eq. 1 ) then 
                       aut_= (1. - 0.1*abs(pw_))**5
                    else
                       aut_= 1 - 0.5*abs(pw_)
                    endif
                    au_=amax1(aa_,aut_)
                    aw_(i,j,k)=dw_*au_+amax1(fw_,aa_)


                    dn_=(gamma_*xdu_(i)*zdu_(k))/ydnu_(j)
                    vn_=(vc_(i,j,k)*xveu_(i)+vc_(i+1,j,k)*xvwu_(i))/(xvwu_(i)+xveu_(i))
                    fn_=vn_*xdu_(i)*zdu_(k)
                    pn_=fn_/dn_
                    if( ischeme_ .eq. 1 ) then 
                       aut_= (1. - 0.1*abs(pn_))**5
                    else
                       aut_= 1 - 0.5*abs(pn_)
                    endif
                    au_=amax1(aa_,aut_)
                    an_(i,j,k)=dn_*au_+amax1(-fn_,aa_)
                    
                    ds_=(gamma_*xdu_(i)*zdu_(k))/ydsu_(j)
                    vs_=(vc_(i,j-1,k)*xveu_(i)+vc_(i+1,j-1,k)*xvwu_(i))/(xvwu_(i)+xveu_(i))
                    fs_=vs_*xdu_(i)*zdu_(k)
                    ps_=fs_/ds_
                    if( ischeme_ .eq. 1 ) then 
                       aut_= (1. - 0.1*abs(ps_))**5
                    else
                       aut_= 1 - 0.5*abs(ps_)
                    endif
                    au_=amax1(aa_,aut_)
                    as_(i,j,k)=ds_*au_+amax1(fs_,aa_)
                    !c-->

                    dt_=(gamma_*xdu_(i)*ydu_(j))/zdtu_(k)
                    wt_=(wc_(i,j,k)*xweu_(i)+wc_(i+1,j,k)*xwwu_(i))/(xwwu_(i)+xweu_(i))
                    ft_=wt_*xdu_(i)*ydu_(j)
                    pt_=ft_/dt_
                    if( ischeme_ .eq. 1 ) then 
                       aut_= (1. - 0.1*abs(pt_))**5
                    else
                       aut_= 1 - 0.5*abs(pt_)
                    endif
                    au_=amax1(aa_,aut_)
                    at_(i,j,k)=dt_*au_+amax1(-ft_,aa_)

                    db_=(gamma_*xdu_(i)*ydu_(j))/zdbu_(k)
                    wb_=(wc_(i,j,k-1)*xweu_(i)+wc_(i+1,j,k-1)*xwwu_(i))/(xwwu_(i)+xweu_(i))
                    fb_=wb_*xdu_(i)*ydu_(j)
                    pb_=fb_/db_
                    if( ischeme_ .eq. 1 ) then 
                       aut_= (1. - 0.1*abs(pb_))**5
                    else
                       aut_= 1 - 0.5*abs(pb_)
                    endif
                    au_=amax1(aa_,aut_)
                    ab_(i,j,k)=db_*au_+amax1(fb_,aa_)

                    ap_(i,j,k)=ae_(i,j,k)+aw_(i,j,k)+an_(i,j,k)+as_(i,j,k) + at_(i,j,k) + ab_(i,j,k) + &
                         (xdu_(i)*ydu_(j)*zdu_(k)/dtau_)
                    st1_ = (t_(i,j,k)*xveu_(i) + t_(i+1,j,k)*xvwu_(i))/(xvwu_(i)+xveu_(i))
                    st_ = st1_*(grast_*(prndl_**2))*xdu_(i)*ydu_(j)*zdu_(k)*sin(phir1_)
                    sc1_ = (c_(i,j,k)*xveu_(i) + c_(i+1,j,k)*xvwu_(i))/(xvwu_(i)+xveu_(i))
                    sc_ = sc1_*(grasm_*(prndl_**2))*xdu_(i)*ydu_(j)*zdu_(k)*sign_*sin(phir1_)
                    avs_(i,j,k) = st_ + sc_ 
                    avo_(i,j,k) =( xdu_(i)*ydu_(j)*zdu_(k)/dtau_)*uc_p(i,j,k)
  
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine ucoef1


subroutine vcoef1(m,n,p,uc_,vc_p,vc_,wc_,t_,c_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,prndl_,grasm_,grast_,dtau_, &
     sign_,ischeme_,phir1_,phir2_,phir3_,xdv_,xdev_,xdwv_,ydv_,ydnv_, &
     ydsv_,yusv_,yunv_,ywsv_,ywnv_,zdv_,zdtv_,zdbv_,be_,bw_,bn_,bs_,bt_,bb_,bp_,bvo_,bvs_,my_id)

  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,ischeme_,my_id
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_,vc_,vc_p,wc_,t_,c_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: be_,bw_,bn_,bs_,bt_,bb_,bp_,bvo_,bvs_
  real prndl_,grasm_,grast_,dtau_,phir1_,phir2_,phir3_,de_,ue_,fe_,pe_,aa_,avt_,av_,gamma_
  real dw_,uw_,fw_,pw_,dn_,vn_,fn_,pn_,ds_,vs_,fs_,ps_,dt_,wt_,ft_,pt_,db_,wb_,fb_,pb_,st_,sc_,st1_,sc1_,sign_
  real, dimension(m) :: xdv_,xdev_,xdwv_
  real, dimension(n) :: ydv_,ydnv_,ydsv_,yusv_,yunv_,ywsv_,ywnv_
  real, dimension(p) :: zdv_,zdtv_,zdbv_
  gamma_= prndl_
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           if ((j .gt. 1) .and. (j .lt. ny_)) then
              do i = sx_,ex_
                 if (i .gt. 1) then
                    de_=(gamma_*ydv_(j)*zdv_(k))/xdev_(i)
                    ue_=(uc_(i,j+1,k)*yusv_(j) + uc_(i,j,k)*yunv_(j))/(yunv_(j) + yusv_(j))
                    fe_=ue_*ydv_(j)*zdv_(k)
                    pe_=fe_/de_
                    aa_=0.
                    if( ischeme_ .eq. 1 ) then 
                       avt_= (1. - 0.1*abs(pe_))**5
                    else
                       avt_= 1 - 0.5*abs(pe_)
                    endif
                    av_=amax1(aa_,avt_)
                    be_(i,j,k)=de_*av_+amax1(-fe_,aa_)
                    !c-->
                    dw_=(gamma_*ydv_(j)*zdv_(k))/xdwv_(i)
                    uw_=(uc_(i-1,j+1,k)*yusv_(j) + uc_(i-1,j,k)*yunv_(j))/(yunv_(j) + yusv_(j))
                    fw_=uw_*ydv_(j)*zdv_(k)
                    pw_=fw_/dw_
                    if( ischeme_ .eq. 1 ) then 
                       avt_= (1. - 0.1*abs(pw_))**5
                    else
                       avt_= 1 - 0.5*abs(pw_)
                    endif
                    av_=amax1(aa_,avt_)
                    bw_(i,j,k)=dw_*av_+amax1(fw_,aa_)


                    dn_=(gamma_*xdv_(i)*zdv_(k))/ydnv_(j)
                    vn_=(vc_(i,j,k) + vc_(i,j+1,k))*0.5
                    fn_=vn_*xdv_(i)*zdv_(k)
                    pn_=fn_/dn_
                    if( ischeme_ .eq. 1 ) then 
                       avt_= (1. - 0.1*abs(pn_))**5
                    else
                       avt_= 1 - 0.5*abs(pn_)
                    endif
                    av_=amax1(aa_,avt_)
                    bn_(i,j,k)=dn_*av_+amax1(-fn_,aa_)
                    
                    ds_=(gamma_*xdv_(i)*zdv_(k))/ydsv_(j)
                    vs_=(vc_(i,j,k) + vc_(i,j-1,k))*0.5
                    fs_=vs_*xdv_(i)*zdv_(k)
                    ps_=fs_/ds_
                    if( ischeme_ .eq. 1 ) then 
                       avt_= (1. - 0.1*abs(ps_))**5
                    else
                       avt_= 1 - 0.5*abs(ps_)
                    endif
                    av_=amax1(aa_,avt_)
                    bs_(i,j,k)=ds_*av_+amax1(fs_,aa_)
                    !c-->

                    dt_=(gamma_*xdv_(i)*ydv_(j))/zdtv_(k)
                    wt_=(wc_(i,j,k)*ywnv_(j)+wc_(i,j+1,k)*ywsv_(j))/(ywnv_(j) + ywsv_(j))
                    ft_=wt_*xdv_(i)*ydv_(j)
                    pt_=ft_/dt_
                    if( ischeme_ .eq. 1 ) then 
                       avt_= (1. - 0.1*abs(pt_))**5
                    else
                       avt_= 1 - 0.5*abs(pt_)
                    endif
                    av_=amax1(aa_,avt_)
                    bt_(i,j,k)=dt_*av_+amax1(-ft_,aa_)

                    db_=(gamma_*xdv_(i)*ydv_(j))/zdbv_(k)
                    wb_=(wc_(i,j,k-1)*ywnv_(j)+wc_(i,j+1,k-1)*ywsv_(j))/(ywnv_(j) + ywsv_(j))
                    fb_=wb_*xdv_(i)*ydv_(j)
                    pb_=fb_/db_
                    if( ischeme_ .eq. 1 ) then 
                       avt_= (1. - 0.1*abs(pb_))**5
                    else
                       avt_= 1 - 0.5*abs(pb_)
                    endif
                    av_=amax1(aa_,avt_)
                    bb_(i,j,k)=db_*av_+amax1(fb_,aa_)

                    bp_(i,j,k) = be_(i,j,k)+bw_(i,j,k)+bn_(i,j,k)+bs_(i,j,k) + bt_(i,j,k) + bb_(i,j,k) + &
                         (xdv_(i)*ydv_(j)*zdv_(k)/dtau_)
                    st1_ = (t_(i,j,k)*yunv_(j) + t_(i,j+1,k)*yusv_(j))/(yusv_(j) + yunv_(j))
                    st_ = st1_*(grast_*(prndl_**2))*xdv_(i)*ydv_(j)*zdv_(k)*sin(phir2_)
                    sc1_ = (c_(i,j,k)*yunv_(j) + c_(i,j+1,k)*yusv_(j))/(yusv_(j) + yunv_(j))
                    sc_ = sc1_*(grasm_*(prndl_**2))*xdv_(i)*ydv_(j)*zdv_(k)*sign_*sin(phir2_)
                    bvs_(i,j,k) = st_ + sc_ 
                    bvo_(i,j,k) =( xdv_(i)*ydv_(j)*zdv_(k)/dtau_)*vc_p(i,j,k)
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine vcoef1


subroutine wcoef1(m,n,p,uc_,vc_,wc_p,wc_,t_,c_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,prndl_,grasm_,grast_,dtau_, &
     sign_,ischeme_,phir1_,phir2_,phir3_,xdw_,xdew_,xdww_,ydw_,ydnw_, &
     ydsw_,zdw_,zdtw_,zdbw_,zubw_,zutw_,zvbw_,zvtw_,ce_,cw_,cn_,cs_,ct_,cb_,cp_,cvo_,cvs_,my_id)
    
  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,ischeme_,my_id
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_,vc_,wc_,wc_p,t_,c_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: ce_,cw_,cn_,cs_,ct_,cb_,cp_,cvo_,cvs_
  real prndl_,grasm_,grast_,dtau_,phir1_,phir2_,phir3_,de_,ue_,fe_,pe_,aa_,awt_,aw_,gamma_
  real dw_,uw_,fw_,pw_,dn_,vn_,fn_,pn_,ds_,vs_,fs_,ps_,dt_,wt_,ft_,pt_,db_,wb_,fb_,pb_,st_,sc_,st1_,sc1_,sign_
  real, dimension(m) :: xdw_,xdew_,xdww_
  real, dimension(n) :: ydw_,ydnw_,ydsw_
  real, dimension(p) :: zdw_,zdtw_,zdbw_,zubw_,zutw_,zvbw_,zvtw_
  gamma_= prndl_
  do k = sz_,ez_
     if ((k .gt. 1) .and. (k .lt. nz_)) then
        do j = sy_,ey_
           if (j .gt. 1) then
              do i = sx_,ex_
                 if (i .gt. 1) then
                    de_ = (gamma_*zdw_(k)*ydw_(j))/xdew_(i)
                    ue_ = (uc_(i,j,k)*zutw_(k) + uc_(i,j,k+1)*zubw_(k))/(zutw_(k) + zubw_(k))
                    fe_ = ue_*ydw_(j)*zdw_(k)
                    pe_=fe_/de_
                    aa_=0.
                    if (ischeme_ .eq. 1) then
                       awt_ = (1. - 0.1*abs(pe_))**5
                    else
                       awt_ = 1. - 0.5*abs(pe_)
                    endif
                    aw_=amax1(aa_,awt_)
                    ce_(i,j,k)=de_*aw_+amax1(-fe_,aa_)

                    ! if( islip_ .eq. 1 .and. i .eq. nx ) then
                    !    ce(i,j,k) = 0
                    ! endif

                    !c--   >
                    dw_ = (gamma_*zdw_(k)*ydw_(j))/xdww_(i)
                    uw_ = (uc_(i-1,j,k)*zutw_(k) + uc_(i-1,j,k+1)*zubw_(k))/(zutw_(k) + zubw_(k))
                    fw_ = uw_*ydw_(j)*zdw_(k)
                    pw_=fw_/dw_
                    aa_=0.
                    if (ischeme_ .eq. 1) then
                       awt_ = (1. - 0.1*abs(pw_))**5
                    else
                       awt_ = 1. - 0.5*abs(pw_)
                    endif
                    aw_=amax1(aa_,awt_)
                    cw_(i,j,k)=dw_*aw_+amax1(fw_,aa_)

                    ! if( islip .eq. 1 .and. i .eq. 2 ) then
                    !    cw(i,j,k) = 0
                    ! endif
                    ! c--   >
                    dn_ = (gamma_*xdw_(i)*zdw_(k))/ydnw_(j)
                    vn_ = (vc_(i,j,k)*zvtw_(k) + vc_(i,j,k+1)*zvbw_(k))/(zvtw_(k) + zvbw_(k))
                    fn_ = vn_*xdw_(i)*zdw_(k)
                    pn_ = fn_/dn_
                    if (ischeme_ .eq. 1) then
                       awt_ = (1. - 0.1*abs(pn_))**5
                    else
                       awt_ = 1. - 0.5*abs(pn_)
                    endif
                    aw_=amax1(aa_,awt_)
                    cn_(i,j,k)=dn_*aw_+amax1(-fn_,aa_)

                    ! if( islip .eq. 1 .and. j .eq. ny ) then
                    !    cn(i,j,k) = 0
                    ! endif
                    !c---  > 
                    ds_ = (gamma_*xdw_(i)*zdw_(k))/ydsw_(j)
                    vs_ = (vc_(i,j-1,k)*zvtw_(k) + vc_(i,j-1,k+1)*zvbw_(k))/(zvtw_(k) + zvbw_(k))
                    fs_ = vs_*xdw_(i)*zdw_(k)
                    ps_ = fs_/ds_
                    if (ischeme_ .eq. 1) then
                       awt_ = (1. - 0.1*abs(ps_))**5
                    else
                       awt_ = 1. - 0.5*abs(ps_)
                    endif
                    aw_=amax1(aa_,awt_)
                    cs_(i,j,k)=ds_*aw_+amax1(fs_,aa_)

                    ! if( islip .eq. 1 .and. j .eq. 2 ) then
                    !    cs(i,j,k) = 0
                    ! endif
                    ! c---- >
                    dt_ = (gamma_*xdw_(i)*ydw_(j))/zdtw_(k)
                    wt_ = (wc_(i,j,k) + wc_(i,j,k+1))*0.5
                    ft_ = wt_*xdw_(i)*ydw_(j)
                    pt_ = ft_/dt_
                    if (ischeme_ .eq. 1) then
                       awt_ = (1. - 0.1*abs(pt_))**5
                    else
                       awt_ = 1. - 0.5*abs(pt_)
                    endif
                    aw_=amax1(aa_,awt_)
                    ct_(i,j,k)=dt_*aw_+amax1(-ft_,aa_)

                    db_ = (gamma_*xdw_(i)*ydw_(j))/zdbw_(k)
                    wb_ = (wc_(i,j,k) + wc_(i,j,k-1))*0.5
                    fb_ = wb_*xdw_(i)*ydw_(j)
                    pb_ = fb_/db_
                    if (ischeme_ .eq. 1) then
                       awt_ = (1. - 0.1*abs(pb_))**5
                    else
                       awt_ = 1. - 0.5*abs(pb_)
                    endif
                    aw_ = amax1(aa_,awt_)
                    cb_(i,j,k) = db_*aw_ + amax1(fb_,aa_)

                    cp_(i,j,k) = ce_(i,j,k) + cw_(i,j,k) + cn_(i,j,k) + cs_(i,j,k)+ &
                         ct_(i,j,k) + cb_(i,j,k) + (xdw_(i)*ydw_(j)*zdw_(k)/dtau_)
                    st1_ = (t_(i,j,k)*zutw_(k) + t_(i,j,k+1)*zubw_(k))/(zutw_(k) + zubw_(k))
                    st_ = st1_*(grast_*(prndl_**2))*xdw_(i)*ydw_(j)*zdw_(k)*cos(phir3_)
                    sc1_ = (c_(i,j,k)*zutw_(k)+c_(i,j,k+1)*zubw_(k))/(zutw_(k)+zubw_(k))
                    sc_ = sc1_*(grasm_*(prndl_**2))*xdw_(i)*ydw_(j)*zdw_(k)*sign_*cos(phir3_)
                    cvs_(i,j,k) = st_ + sc_ 
                    cvo_(i,j,k) = (xdw_(i)*ydw_(j)*zdw_(k)/dtau_ )*wc_p(i,j,k)
!                    print *, i,j,k,cvs_(i,j,k),cvo_(i,j,k),zutw_(k),zubw_(k),st1_,sc1_,my_id
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine wcoef1
