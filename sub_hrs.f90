subroutine sub_hrs(lh, hmin, hmax, c, uc, indh, h)
    
    use mod_gvar
    use omp_lib
    use mod_functions
    
    implicit none   
    
    integer, intent(in) :: lh
    
    double precision, intent (in) :: hmin, hmax, c, uc 
    
    integer, intent(out) :: indh
    
    double precision, intent(out) :: h
    
    integer :: ih
    
    double precision :: hgrid(lh), uh_vec(lh), xx, yy, taxss, taxw, mrgtaxw, mrgtaxss, intra_RHS(lh), diff(lh)
    
        ! hgrid
        hgrid(1) = hmin
        hgrid(lh) = hmax
        do ih = 2,lh-1
            hgrid(ih) = hgrid(1) + (hgrid(lh) - hgrid(1)) / (lh-1) * (ih-1)
        end do
                    
        do ih = 1,lh
            call sub_util(c,hgrid(ih),xx,xx,uh_vec(ih))
            yy = w * eprof_age_ablt(age,n) * hgrid(ih)
            taxss = tau_ss * min(y_bar, yy)
            call sub_tax_hsv(prog, yy-0.5*taxss,taxw,mrgtaxw)
            mrgtaxss = 0.
            mrgtaxss = tau_ss * (yy > y_bar)
            intra_RHS(ih) = uh_vec(ih) / (w * eprof_age_ablt(age,n) * eps(is,age,n) * (1.- mrgtaxw - mrgtaxss))
            diff(ih) = uc / (1-tau_c) - intra_RHS(ih)
        end do
        indh = minloc(abs(diff), dim = 1)
                    
        h = hgrid(indh) 

end subroutine

   
    