module mod_functions
    
    use mod_gvar
    
    contains
    
    !----------------------------------------------------------------------------------------------
    subroutine sub_util(c,h,u,uc,ul)
    
    implicit none
        
        double precision, INTENT(in)  :: c,h
        double precision, INTENT(OUT) :: u, uc, ul 
        
        u  = c**(1.-sigma) / (1.0-sigma) + B * (1.-h)**(1.+1./xi) / (1.+1./xi)
        uc = c**(-sigma) 
        ul = B * (1.-h)**(1./xi) 
        
    end subroutine
    
    !----------------------------------------------------------------------------------------------
    subroutine sub_util_inverse(deriv_V, h, c)
    
    implicit none
        
        double precision, INTENT(in)  :: deriv_V, h
        double precision, INTENT(OUT) :: c 
        
        c = (deriv_V)**(-1./sigma)
        
    end subroutine
    
    !----------------------------------------------------------------------------------------------
    subroutine sub_weight(lg, grid, x, ind, wgt)
    
    implicit none        
    
    integer, intent(in) :: lg
    double precision, dimension(:), intent(in) :: grid(lg)
    double precision, intent(in) :: x
    
    integer :: ind
    
    double precision :: wgt    
    
        ! x is within grid(1) and grid(lg)
        x = max(x, grid(1)  + 1e-8)
        x = min(x, grid(lg) - 1e-8)
        
        do ix = 1,lg    
            diff(ix) = abs(x - grid(ix))
        end do

        iix = minloc(diff, dim = 1)

        if ( x .ge. grid(iix) ) then
            ind = iix
        elseif ( x < grid(iix) ) then
            ind = iix-1
        end if
        
        wgt = (grid(inda+1) - x) / (grid(inda+1) - grid(inda))
    
    end subroutine
    
    !----------------------------------------------------------------------------------------------
    
    subroutine sub_tax_hsv(prog, y, tax, mrgtax)
    
        implicit none

        double precision, INTENT(in)  :: prog, y
        double precision, INTENT(OUT) :: y, tax, mrgtax
            
        if (prog == 0) then
            ! flat tax
            tax = tau_w * y
            mrgtax = tau_w
        else
            ! progressive tax
            tax = y - lambda_l * y**(1.- tau_w) 
            mrgtax = 1 - lambda_l * (1.- tau_w) * y**(-tau_w)
        end if
        
    end subroutine   ! sub_hsvtax
    
    
    !!----------------------------------------------------------------------------------------------
    !subroutine sub_init_dist
    !    
    !    implicit none
    !    
    !    double precision :: pass_init_low(822,3), pass_init_high(372,3), diff_a(la)
    !    integer :: inda, iia, ja
    !    
    !    
    !    ! initial distribution over assets, prod & skill
    !    !OPEN(3,file=fileplace2//"/initial_asset0.txt",status='old',action='READ') 
    !    OPEN(3,file=fileplace5//"/initial_asset0.txt",status='old',action='READ') 
    !    do ia = 1,822
    !        !         shock is             assests ia           probability
    !        READ(3,*) pass_init_low(ia,1), pass_init_low(ia,2), pass_init_low(ia,3) 
    !    end do
    !
    !    
    !    !OPEN(3,file=fileplace2//"/initial_asset1.txt",status='old',action='READ') 
    !    OPEN(3,file=fileplace5//"/initial_asset1.txt",status='old',action='READ') 
    !    do ia = 1,372
    !        READ(3,*) pass_init_high(ia,1), pass_init_high(ia,2), pass_init_high(ia,3)
    !    end do
    !
    !    
    !    ! annual income 60k, 2 year GDP 120k
    !    pass_init_low( :,2) = pass_init_low( :,2) / (iyear *  60000.0) * Y0
    !    pass_init_high(:,2) = pass_init_high(:,2) / (iyear *  60000.0) * Y0
    !
    !
    !    prob_pol = 0.0   
    !    
    !    n = 1
    !    age = 1
    !    ih = ih0_low
    !    
    !    do ia = 1,822
    !        do is = 1,ls            
    !            
    !            if (abs(pass_init_low(ia,1) - is) < 1e-4) then ! the following distriubtion over asset belongs to shock s
    !                
    !                ! asset position, pass_init_low(i,2)
    !                call sub_weight(agrid, pass_init_low(ia,2), ind, wgt)
    !                
    !                ! starting with bequest capital
    !                prob_pol(ind,  is, ih, age, n) = pass_init_low(ia,3) *      wgt  + prob_pol(ind,  is, ih, age, n)
    !                prob_pol(ind+1,is, ih, age, n) = pass_init_low(ia,3) * (1.0-wgt) + prob_pol(ind+1,is, ih, age, n)
    !                
    !                prob_pol(ind   is, ih, age, n) = max(0.0, prob_pol(ind,  is, ih, age, n))
    !                prob_pol(ind+1,is, ih, age, n) = max(0.0, prob_pol(ind+1,is, ih, age, n))
    !                
    !            end if
    !        end do
    !    end do
    !        
    !    do is = 1,ls
    !        prob_pol(1:la, is, ih, age, n) = prob_pol(1:la, is, ih, age, n) * pstar(is,1)
    !    end do 
    !    prob_pol(1:la, 1:ls, ih, age, 1) = prob_pol(1:la, 1:ls, ih, age, 1) / sum(prob_pol(1:la, 1:ls, ih, age, 1))
    !
    !
    !    n = 2
    !    age = 1
    !    ih = ih0_high
    !    
    !    do ia = 1,372
    !        do is = 1,ls
    !            
    !            if (pass_init_high(ia,1) == is) then ! the following distriubtion over asset belongs to shock s!                    
    !
    !                ! asset position, pass_init_high(i,2)               
    !                call sub_weight(agrid, pass_init_low(ia,2), ind, wgt)
    !                                    
    !                ! starting with bequest capital
    !                prob_pol(ind,  is, ih, age, n) = pass_init_high(ia,3) *      wgt  + prob_pol(ind,  is, ih, age, n)
    !                prob_pol(ind+1,is, ih, age, n) = pass_init_high(ia,3) * (1.0-wgt) + prob_pol(ind+1,is, ih, age, n)
    !                
    !                prob_pol(ind,  is, ih, age, n) = max(0.0, prob_pol(ind,  is, ih, age, n))
    !                prob_pol(ind+1,is, ih, age, n) = max(0.0, prob_pol(ind+1,is, ih, age, n))
    !                
    !            end if
    !        end do
    !    end do
    !
    !    do is = 1,ls
    !        prob_pol(1:la, is, ih, age, n) = prob_pol(1:la, is, ih, age, n) * pstar(is,1)
    !    end do         
    !    prob_pol(1:la, 1:ls, ih, age, 2) = prob_pol(1:la, 1:ls, ih, age, 2) / sum(prob_pol(1:la, 1:ls, ih, age, 2))
    !    
    !
    !end subroutine

        
end module