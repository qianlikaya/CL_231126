subroutine sub_con_hrs(lh, hmin, hmax, ap, EVhata, indh, h)
    
    use mod_gvar
    use omp_lib
    use mod_functions
    
    implicit none   
    
    integer, intent(in) :: lh
    
    double precision, intent (in) :: hmin, hmax, ap, wgta, Vhata, Vhata_1
    
    integer, intent(out) :: indh
    
    double precision, intent(out) :: h
    
    integer :: ih
    
    double precision :: hgrid(lh), xx, yy, taxss, taxw, c, V(lh)
    
        ! hgrid
        hgrid(1) = hmin
        hgrid(lh) = hmax
        do ih = 2,lh-1
            hgrid(ih) = hgrid(1) + (hgrid(lh) - hgrid(1)) / (lh-1) * (ih-1)
        end do
                    
        do ih = 1,lh
            yy = w * eprof_age_ablt(age,n) * eps(is,age,n) * hgrid(ih)
            taxss = tau_ss * min(y_bar, yy)
            call sub_tax_hsv(prog, yy-0.5*taxss,taxw, xx)
                        
            c = ( (1+(1-tau_a)*r) * (agrid(ia) + bequest) + yy - taxw - taxss + lump - ap ) / (1+tau_c) 
                        
            call sub_util(c, hgrid(ih), u, xx, xx)
                        
            V(ih) = u + EVhata
                        
        end do
                    
        indh = maxloc(V, dim=1)
        h = hgrid(indh)

end subroutine

   
    