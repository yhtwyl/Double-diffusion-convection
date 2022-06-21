subroutine prcoff(m,n,p,dum1_,dum2_,dum3_,ucp_,vcp_,wcp_,pe_,pw_,pn_,ps_,pt_,pbtm_,pp_,pb_,xdt_,ydt_,zdt_,nx_,ny_,nz_,sx_,ex_,sy_,&
     ey_,sz_,ez_,my_id)
!  parameter(m=22,n=22)

  integer m,n,p,i,j,k,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_,my_id
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: dum1_,dum2_,dum3_,ucp_,vcp_,wcp_,pe_,pw_,pn_,ps_,pt_,pbtm_,pp_,pb_
  real, dimension(m) :: xdt_
  real, dimension(n) :: ydt_
  real, dimension(p) :: zdt_
  
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           if (j .gt. 1) then
              do i = sx_,ex_
                 if (i .gt. 1) then
                    if ((i .ne. 2 ) .and. (i .ne. nx_)) then
                       pe_(i,j,k)=(ydt_(j)*ydt_(j)*zdt_(k)*zdt_(k))/dum1_(i,j,k)
                       pw_(i,j,k)=(ydt_(j)*ydt_(j)*zdt_(k)*zdt_(k))/dum1_(i-1,j,k)
                    elseif (i .eq. 2) then
                       pe_(i,j,k)=(ydt_(j)*ydt_(j)*zdt_(k)*zdt_(k))/dum1_(i,j,k)
                       pw_(i,j,k)=0.
                    elseif (i .eq. nx_) then
                       pe_(i,j,k)=0.
                       pw_(i,j,k)=(ydt_(j)*ydt_(j)*zdt_(k)*zdt_(k))/dum1_(i-1,j,k)
                    endif

                    if ((j .ne. 2) .and. (j .ne. ny_)) then
                       pn_(i,j,k)=(xdt_(i)*xdt_(i)*zdt_(k)*zdt_(k))/dum2_(i,j,k)
                       ps_(i,j,k)=(xdt_(i)*xdt_(i)*zdt_(k)*zdt_(k))/dum2_(i,j-1,k)
                    elseif (j .eq. 2 ) then
                       pn_(i,j,k)=(xdt_(i)*xdt_(i)*zdt_(k)*zdt_(k))/dum2_(i,j,k)
                       ps_(i,j,k)=0.
                    else if(j .eq. ny_ ) then
                       pn_(i,j,k)=0.
                       ps_(i,j,k)=(xdt_(i)*xdt_(i)*zdt_(k)*zdt_(k))/dum2_(i,j-1,k)
                    endif

                    if ((k .ne. 2) .and. (k .ne. nz_)) then
                       pt_(i,j,k)=(xdt_(i)*xdt_(i)*ydt_(j)*ydt_(j))/dum3_(i,j,k)
                       pbtm_(i,j,k)=(xdt_(i)*xdt_(i)*ydt_(j)*ydt_(j))/dum3_(i,j,k-1)
                    elseif (k .eq. 2) then
                       pt_(i,j,k)=(xdt_(i)*xdt_(i)*ydt_(j)*ydt_(j))/dum3_(i,j,k)
                       pbtm_(i,j,k)=0.
                    else if(k .eq. nz_) then
                       pt_(i,j,k)=0.
                       pbtm_(i,j,k)=(xdt_(i)*xdt_(i)*ydt_(j)*ydt_(j))/dum3_(i,j,k-1)
                    endif
                    pp_(i,j,k)= pe_(i,j,k) + pw_(i,j,k) + pn_(i,j,k) + ps_(i,j,k) + pt_(i,j,k) + pbtm_(i,j,k)
                    pb_(i,j,k)=(ucp_(i-1,j,k)-ucp_(i,j,k))*ydt_(j)*zdt_(k)+ (vcp_(i,j-1,k)-vcp_(i,j,k))*xdt_(i)*zdt_(k) + & 
                         (wcp_(i,j,k-1)-wcp_(i,j,k))*xdt_(i)*ydt_(j)
!                    print *,i,j,k, pe_(i,j,k), pw_(i,j,k), pn_(i,j,k), ps_(i,j,k), pt_(i,j,k), pbtm_(i,j,k),pp_(i,j,k), &
!                         pb_(i,j,k), my_id
                 endif
              enddo
           endif
        enddo
     endif
  enddo
end subroutine prcoff

subroutine presrx(mp,prs_,pe_,pw_,pn_,ps_,pt_,pbtm_,pp_,pb_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_)
  !  parameter(m=22)
  integer mp,i,i_,j,k,nx_,ny_,nz_,start_,sx_,ex_,sy_,ey_,sz_,ez_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: prs_,pe_,pw_,pn_,ps_,pt_,pbtm_,pp_,pb_
  real, dimension(mp) :: a__,b__,c__,d__,x__
  start_ = sx_ - 2;
  do k = sz_,ez_
     if (k .gt. 1) then
        do j = sy_,ey_
           i_ = 0;
           if (j .gt. 1) then
              do i=sx_-1,ex_+1
                 if ((i .gt. 1) .and. (i .lt. nx_+1)) then
                    i_ = i_ + 1;
                    d__(i_) = pn_(i,j,k)*prs_(i,j+1,k) + ps_(i,j,k)*prs_(i,j-1,k) + pt_(i,j,k)*prs_(i,j,k+1) + &
                         pbtm_(i,j,k)*prs_(i,j,k-1) + pb_(i,j,k);
                    b__(i_) = pp_(i,j,k);
                    a__(i_) = -pw_(i,j,k)
                    c__(i_) = -pe_(i,j,k)

                    if (i .eq. 2) then
                       a__(i_) = 0.0;
                       start_ = 1;
                    endif
                    if (i .eq. nx_) c__(i_) = 0.0;
                    if ((i .eq. sx_-1) .or. (i .eq. ex_+1)) then
                       a__(i_) = 0.0;c__(i_)= 0.0;b__(i_) = 1.0;d__(i_)=prs_(i,j,k);
                    endif
                 endif
              enddo
              call tridag(a__,b__,c__,d__,i_,x__,mp,'presrx')
              do i=1,i_
                 prs_(start_ + i,j,k) = x__(i);
                 !           print *,'creator prs_',start_+i,j,prs_(start_+i,j)
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine presrx

subroutine presry(mp,prs_,pe_,pw_,pn_,ps_,pt_,pbtm_,pp_,pb_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_)
  integer mp,i,j_,j,k,nx_,ny_,nz_,start_,end_,sx_,ex_,sy_,ey_,sz_,ez_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: prs_,pe_,pw_,pn_,ps_,pt_,pbtm_,pp_,pb_
  real, dimension(mp) :: a__,b__,c__,d__,x__ !! serious issue for aspect ratio, either make it equal to larger one
  start_ = sy_ - 2;
  do k = sz_,ez_
     if (k .gt. 1) then
        do i = sx_,ex_
           j_ = 0;
           if (i .gt. 1) then
              do j=sy_-1,ey_+1
                 if ((j .gt. 1) .and. (j .lt. ny_+1)) then
                    j_ = j_+1;
                    d__(j_) = pe_(i,j,k)*prs_(i+1,j,k) + pw_(i,j,k)*prs_(i-1,j,k) + pt_(i,j,k)*prs_(i,j,k+1) + &
                         pbtm_(i,j,k)*prs_(i,j,k-1) + pb_(i,j,k);
                    b__(j_) = pp_(i,j,k);
                    a__(j_) = -ps_(i,j,k)
                    c__(j_) = -pn_(i,j,k)
                    
                    if (j .eq. 2) then
                       a__(j_) = 0.0;
                       start_ = 1;
                    endif
                    if (j .eq. ny_) then
                       c__(j_) = 0.0; 
                    endif
                    if ((j .eq. sy_-1) .or. (j .eq. ey_+1)) then
                       a__(j_) = 0.0;c__(j_)= 0.0;b__(j_) = 1.0;d__(j_)=prs_(i,j,k);
                    endif
                 endif
              enddo
              call tridag(a__,b__,c__,d__,j_,x__,mp,'presry')
              do j=1,j_
                 prs_(i,start_ + j,k) = x__(j);
                 !           print *,'creator prs_',start_+i,j,prs_(start_+i,j)
              enddo
           endif
        enddo
     endif
  enddo
  return
end subroutine presry

subroutine presrz(mp,prs_,pe_,pw_,pn_,ps_,pt_,pbtm_,pp_,pb_,nx_,ny_,nz_,sx_,ex_,sy_,ey_,sz_,ez_)
!  parameter(m=22)
  integer mp,i,j,k_,k,nx_,ny_,nz_,start_,end_,sx_,ex_,sy_,ey_,sz_,ez_
  real, dimension(sx_-1:ex_+1,sy_-1:ey_+1,sz_-1:ez_+1) :: prs_,pe_,pw_,pn_,ps_,pt_,pbtm_,pp_,pb_
  real, dimension(mp) :: a__,b__,c__,d__,x__ !! serious issue for aspect ratio, either make it equal to larger one
  start_ = sz_ - 2;
  do j = sy_,ey_
     if (j .gt. 1) then
        do i = sx_,ex_
           k_ = 0;
           if (i .gt. 1) then
              do k = sz_-1,ez_+1
                 if ((k .gt. 1) .and. (k .lt. nz_+1)) then
                    k_ = k_+1;
                    d__(k_) = pe_(i,j,k)*prs_(i+1,j,k) + pw_(i,j,k)*prs_(i-1,j,k) + pn_(i,j,k)*prs_(i,j+1,k) + &
                         ps_(i,j,k)*prs_(i,j-1,k) + pb_(i,j,k);
                    b__(k_) = pp_(i,j,k);
                    a__(k_) = -pbtm_(i,j,k)
                    c__(k_) = -pt_(i,j,k)
                    
                    if (k .eq. 2) then
                       a__(k_) = 0.0;
                       start_ = 1;
                    endif
                    if (k .eq. nz_) then
                       c__(k_) = 0.0; 
                    endif
                    if ((k .eq. sz_-1) .or. (k .eq. ez_+1)) then
                       a__(k_) = 0.0;c__(k_)= 0.0;b__(k_) = 1.0;d__(k_)=prs_(i,j,k);
                    endif
                 endif
              enddo
              call tridag(a__,b__,c__,d__,k_,x__,mp,'presrz')
              do k = 1,k_
                 prs_(i,j,start_+k) = x__(k);
!           print *,'creator prs_',i,j,start_+j,prs_(i,start_+j)
              enddo
           endif
        enddo
     endif
  enddo
end subroutine presrz

