subroutine screen(grast_,grasm_,tau_,dtau_,iter_,itr_,itau_,rmmass_,nxyz_,enbal_,conbal_,umbal_,vmbal_,wmbal_)
  real grast_,grasm_,tau_,dtau_,rmmass_,enbal_,conbal_,umbal_,vmbal_,wmbal_
  integer iter_,itr_,itau_,nxyz_
  write(*,*)'---------------------------------------------'
  write(*,*)
  write(*,*)'grast =  ',grast_,'   grasm =  ',grasm_
  write(*,*)'iter =  ',iter_,'   rmmass =   ',rmmass_
  write(*,*)'nxyz = ',nxyz_,'    enbal =     ',enbal_
  write(*,*)'conbal   = ',conbal_,'  umbal = ',umbal_
  write(*,*)' vmbal =    ',vmbal_ ,' wmbal = ',wmbal_
  write(*,*)' tau  =   ',tau_
  write(*,*)'itr    =',itr_,'    itau =   ',itau_
  write(*,*)'--------------------------------------------'  
end subroutine screen
