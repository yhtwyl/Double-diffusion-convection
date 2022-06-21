subroutine species_coeff(m,n,p,uc_,vc_,wc_,t_p,aet_,awt_,ant_,ast_,atop_,abtm_,apt_,apot_,xdet_,xdt_,xdwt_,ydt_,ydnt_,ydst_&
     ,zdt_,zdbt_,zdtt_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,prndl_,schmdt_,dtau_,ischeme_,ibgrad_,itgrad_,who)

  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,ischeme_,ibgrad_,itgrad_
  real, dimensIon(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: uc_,vc_,wc_,t_p
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: aet_,awt_,ant_,ast_,atop_,abtm_,apt_,apot_
  real dtau_,de_,fe_,pe_,aa_,att_,at_,schmdt_,prndl_,gamma_
  real dw_,fw_,pw_,dn_,fn_,pn_,ds_,fs_,ps_,dt_,ft_,pt_,db_,fb_,pb_
  real, dimension(m) :: xdet_,xdt_,xdwt_
  real, dimension(n) :: ydt_,ydnt_,ydst_
  real, dimension(p) :: zdt_,zdbt_,zdtt_
  character(2) who
  gamma_ = 1.0
  if (who .eq. 'c_') gamma_ = prndl_/schmdt_;

  do k = sz_,ez_
     if (k .gt. 1) then
        do j=sy_,ey_
           if (j .gt. 1) then
              do i = sx_,ex_
                 if (i .gt. 1) then
                    de_=(gamma_*ydt_(j)*zdt_(k))/xdet_(i)
                    fe_=uc_(i,j,k)*ydt_(j)*zdt_(k)
                    pee_=fe_/de_
                    aa_ = 0
                    if( ischeme_ .eq. 1 ) then
                       att_ = ( 1. - 0.1*abs(pee_))**5
                    else
                       att_ = 1. - 0.5*abs(pee_)
                    endif
                    at_=amax1(aa_,att_)
                    aet_(i,j,k)=de_*at_+amax1(-fe_,aa_)
                    if( i .eq. nx_ ) then !! WE ASSUME HERE THAT SIDE WALLS ARE AT ZERO GRADIENT OF SCALARS !!
                       aet_(i,j,k) = 0.
                    endif
!!!!!!!!!!!!!!!!!!c--   >

                    dw_=(gamma_*ydt_(j)*zdt_(k))/xdwt_(i)
                    fw_=uc_(i-1,j,k)*ydt_(j)*zdt_(k)
                    pww_=fw_/dw_
                    if( ischeme_ .eq. 1 ) then
                       att_ = ( 1. - 0.1*abs(pww_))**5
                    else
                       att_ = 1. - 0.5*abs(pww_)
                    endif
                    at_=amax1(aa_,att_)
                    awt_(i,j,k)=dw_*at_+amax1(fw_,aa_)
                    if(i .eq. 2) then
                       awt_(i,j,k) = 0.
                    endif
!!!!!!!!!!!!!!!!!c--   >

                    dn_=(gamma_*xdt_(i)*zdt_(k))/ydnt_(j)
                    fn_=vc_(i,j,k)*xdt_(i)*zdt_(k)
                    pnn_=fn_/dn_
                    if( ischeme_ .eq. 1 ) then
                       att_ = ( 1. - 0.1*abs(pnn_))**5
                    else
                       att_ = 1. - 0.5*abs(pnn_)
                    endif
                    at_=amax1(aa_,att_)
                    ant_(i,j,k)=dn_*at_+amax1(-fn_,aa_)
                    if( j .eq. ny_ ) then
                       ant_(i,j,k) = 0.
                    endif
!!!!!!!!!!!!!!!!!c--   >
                    
                    ds_=(gamma_*xdt_(i)*zdt_(k))/ydst_(j)
                    fs_=vc_(i,j-1,k)*xdt_(i)*zdt_(k)
                    pss_=fs_/ds_
                    if( ischeme_ .eq. 1 ) then
                       att_ = ( 1. - 0.1*abs(pss_))**5
                    else
                       att_ = 1. - 0.5*abs(pss_)
                    endif
                    at_=amax1(aa_,att_)
                    ast_(i,j,k)=ds_*at_+amax1(fs_,aa_)
                    if( j .eq. 2 ) then
                       ast_(i,j,k) = 0.0
                    endif
                    
!!!!!!!!!!!!!!!!!c--   >
                    
                    dt_=(gamma_*xdt_(i)*ydt_(j))/zdtt_(k)
                    ft_=wc_(i,j,k)*xdt_(i)*ydt_(j)
                    ptt_=ft_/dt_
                    if( ischeme_ .eq. 1 ) then
                       att_ = ( 1. - 0.1*abs(ptt_))**5
                    else
                       att_ = 1. - 0.5*abs(ptt_)
                    endif
                    at_=amax1(aa_,att_)
                    atop_(i,j,k)=dt_*at_+amax1(-ft_,aa_)
                    if( k .eq. nz_ ) then
                       atop_(i,j,k) = 0.0
                    endif
                    
!!!!!!!!!!!!!!!!!!c--   >
                    
                    db_=(gamma_*xdt_(i)*ydt_(j))/zdbt_(k)
                    fb_=wc_(i,j,k-1)*xdt_(i)*ydt_(j)
                    pbb_=fb_/db_
                    if( ischeme_ .eq. 1 ) then
                       att_ = ( 1. - 0.1*abs(pbb_))**5
                    else
                       att_ = 1. - 0.5*abs(pbb_)
                    endif
                    at_=amax1(aa_,att_)
                    abtm_(i,j,k)=db_*at_ + amax1(fb_,aa_)
                    if( k .eq. 2 ) then
                       abtm_(i,j,k) = 0.0
                    endif
                    
                    apt_(i,j,k) = aet_(i,j,k) + awt_(i,j,k) + ant_(i,j,k) + ast_(i,j,k) + atop_(i,j,k) + abtm_(i,j,k) +&
                         (xdt_(i)*ydt_(j)*zdt_(k))/dtau_
                    apot_(i,j,k) =(xdt_(i)*ydt_(j)*zdt_(k)/dtau_)*t_p(i,j,k)
                    
!                    print *, i,j,k,aet_(i,j,k),awt_(i,j,k),ant_(i,j,k),ast_(i,j,k),atop_(i,j,k),abtm_(i,j,k), apot_(i,j,k)
                 endif
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine species_coeff
