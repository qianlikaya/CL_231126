subroutine sub_aggregation(Na,Ns, Nt, Nn, x_pol, p_pol, x_mod)
    
    use mod_gvar
    use omp_lib
    
    implicit none   
    
    integer, intent(in) :: Na,Ns, Nt, Nn
    
    double precision, dimension(:), intent (in) :: x_pol(Na,Ns, Nt, Nn), p_pol(Na,Ns, Nt, Nn)
    
    double precision, intent(out) :: x_mod
    
    x_mod = sum(x_pol * p_pol) / sum(p_pol)

end subroutine

   
    