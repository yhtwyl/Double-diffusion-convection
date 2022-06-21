subroutine pressure_calc_first(m_d,s_1,s_2,mp,prs_,a_,b_,c_,d_,e_,f_,g_,h_,na_,nb_,nc_,sa_,ea_,sb_,eb_,sc_,ec_,sx_,ex_,sy_,ey_,sz_,ez_)
  integer m_d,s_1,s_2,mp,na_,nb_,nc_,sa_,ea_,sb_,eb_,sc_,ec_,sx_,ex_,sy_,ey_,sz_,ez_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: prs_,pe_,pw_,pn_,ps_,pt_,pbtm_,pp_,pb_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  start_ = sa_ - 2;
  do s_2 = sc_,ec_
     if (s_2 .gt. 1) then
        do s_1 = sb_,eb_
           m_d_ = 0;
           if (s_1 .gt. 1) then
              do m_d = sa_,ea_
                 if ((m_d .gt. 1) .and. (m_d .lt. na_ + 1)) then
                    m_d_ = m_d_ + 1;
                    d__(m_d_) = pn_(,j,k)*prs_(i,j+1,k) + ps_(i,j,k)*prs_(i,j-1,k) + pt_(i,j,k)*prs_(i,j,k+1) + pbtm_(i,j,k)*prs_(i,j,k-1) + pb_(i,j,k);
                    b__(i_) = pp_(i,j,k);
                    a__(i_) = -pw_(i,j,k)
                    c__(i_) = -pe_(i,j,k)
                    
end subroutine pressure_calc_first
