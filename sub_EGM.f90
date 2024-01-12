subroutine sub_egm(Vhat, Vchat, &
                &  a_pol_egm, c_pol_egm, h_pol_egm, &
                &  V_pol_egm, Vc_pol_egm, &
                &  taxc_pol_egm, taxw_pol_egm, taxa_pol_egm, taxss_pol_egm, sos_pol_egm ) 
    
    !use mkl_service
    use mod_gvar
    use mod_functions
    
    implicit none
    
    ! parameters in / out
    double precision, intent(in)  :: Vhat(la), Vchat(la)
    double precision, intent(out) ::   a_pol_egm(la),    c_pol_egm(la),    h_pol_egm(la), &
                                &      V_pol_egm(la),   Vc_pol_egm(la), &
                                &   taxc_pol_egm(la), taxw_pol_egm(la), taxa_pol_egm(la), taxss_pol_egm(la), sos_pol_egm(la)
    
    double precision :: deriv_Vhat(la), &
    
                        ! cash on hand 
                    &   coh(la), c_temp, yy_temp, sos_temp, taxw_temp, mrgtaxw_temp, taxss_temp, a0(la), &
        
                        ! hours worked
                    &   uc_temp, hmin, hmax, h_vec(la), h_temp, &
                    
                        ! concave / non-concave region
                    &   diff_dv(la), noncv_deriv_Vmax, noncv_deriv_Vmin, &
                    
                        ! global solution
                    &   coh_temp, u_temp, Vhat_temp(la), max_temp, resource, resource1, resource2, & 
                    
                        ! interpolation
                    &   apgrid_pack(la), a0_pack(la), Vhat_pack(la), Vchat_pack(la), a0_cutoff, &
                    &   apgrid_pack_temp(la), a0_pack_temp(la), Vhat_pack_temp(la), Vchat_pack_temp(la), &
                    &   c_glob_first, u_glob_first, &
                    &   c_noncv_ind_min, u_noncv_ind_min, uc_noncv_ind_min, V_noncv_ind_min, Vc_noncv_ind_min, V_glob_first, &
                    &   erro_v, &
                    &   u, uc, un, &
                    &   wgta_vec(la), yy_vec(la), taxw_vec(la), taxss_vec(la), EVhat(la), EVchat(la), &
                    &   xx
                        
    integer :: ih, lh, indh, &
            &   noncv_ind_min, noncv_ind_max, glob_ind(la), ind_jja, glob_ind_first, iter_delta_v, ja, inda_vec(la), jja, la_pack, iter_v, &
            &  noncv_ind_min_gt_0, noncv_ind_max_gt_0
                        
    logical :: diff_mask(la), diff_mask_1(la), mask_glob_ind(la)
    
                     
        ! ===========================================================================================================================================
        ! derivatives 
        ! derivative of Vhat w.r.t. apgrid
        ! the first grid, right derivative
        deriv_Vhat(1) = (Vhat(2) - Vhat(1)) / (apgrid(2) - apgrid(1))   
                            
        ! the middle derivative
        do ja = 2,la-1
            deriv_Vhat(ja) = ( (Vhat(ja+1) - Vhat(ja)) / (apgrid(ja+1) - apgrid(ja)) + (Vhat(ja) - Vhat(ja-1)) / (apgrid(ja) - apgrid(ja-1)) ) / 2.0
        end do
                            
        ! the last point, left derivative
        deriv_Vhat(la) = (Vhat(la) - Vhat(la-1)) / (apgrid(la) - apgrid(la-1))                       
        
        
        ! ===========================================================================================================================================
        ! cash in hand 
        ! Carrol 2006 eq.(14) or Fella 2014 eq.(13) u'(c) = Vhat'(a') = deriv_Vhat(a')
        
        do ja = 1,la
                                
            ! Euler equation
            if ( deriv_Vhat(ja) > 0.0 )  then
                ! seperability, h_temp is unnecesary
                call sub_util_inverse(deriv_Vhat(ja), h_temp, c_temp) ! seperability, h_temp is unnecesary
            else
                c_temp = 1.0e-10
            end if
            
            ! cash in hand
            coh(ja) = (1 + tau_c) * c_temp + apgrid(ja)
            
            ! labor income, pension, tax, transfers
            if (age > tr ) then  ! retired
                
                ! social security benefits and its tax
                yy_temp = replace * earn_ablt(n) 
                call sub_tax_hsv(prog, sos_temp, taxw_temp, xx)
                
                ! social security tax payments
                taxss_temp = 0.
                
            else  
                
                ! hours worked
                ! u'(c) / u'(h) = (1+tau_c) / [partial(earn - taxw - taxss) / partial(h)]
                ! u'(c) = Vhat'(a'), known
                ! partial(earn - taxw - taxss) / partial(h) = w * eprof * eps(is,age,n) - d(taxw)/d(h) - d(taxss)/d(h)
                ! if h \in [0,1], then good
                ! if h < 0, then h = 0.
                ! if h > 1, then h = 1
                uc_temp = deriv_Vhat(ja)
                
                ! discritize hgrid
                ! if ja = 1, meaning hours worked is big, we investigate the entire [0,1] ---> h_vec(ja=1)
                ! for ja > 1, we only investigate a neighborhood of h_vec(ja-1)
                ! if the optimal h is obtained at the first hgrid(1), meaning the region is not wide enough, 
                ! then there might be a jump in the case with fixed labor cost, then we extend the region to include 0
                ! index h_vec for further use
                if (ja == 1) then
                    lh = 100
                    hmin = 0.
                    hmax = 1. - 1e-8
                    call sub_hrs(lh, hmin, hmax, c_temp, uc_temp, indh, h_vec(ja))
                    
                else ! ja > 1
                    ! hgrid
                    lh = 30
                    hmin = max(0., h_vec(ja-1)-0.1)
                    hmax = min(1 - 1e-8, h_vec(ja-1)+0.1)
                    call sub_hrs(lh, hmin, hmax, c_temp, uc_temp, indh, h_vec(ja))
                    
                    ! lower bound not enough, zero hrs might exsit, reduce the lower bound to 0
                    if (indh == 1 .and. hmin > 0) then
                        lh = 30
                        hmin = 0.
                        hmax = h_vec(ja-1)
                        call sub_hrs(lh, hmin, hmax, c_temp, uc_temp, indh, h_vec(ja))
                    end if
                    
                end if
                
                ! earnings, labor tax, social security tax payments
                yy_temp = w * eprof_age_ablt(age,n) * eps(is,age,n) * h_vec(ja)
                taxss_temp  = tau_ss * min(y_bar, yy_temp)
                call sub_tax_hsv(prog, yy_temp-0.5*taxss_temp ,taxw_temp, xx)
                                                                                                      
            end if
                                            
            ! capital of the current period 
            a0(ja) = ( coh(ja) - yy_temp + taxw_temp + taxss_temp - lump ) / (1 + (1-tau_a)*r) - bequest
                        
        end do  

        
        ! ===========================================================================================================================================
        ! concave / non-concave region 
        diff_dv(      2:la) = deriv_Vhat(2:la) - deriv_Vhat(1:la-1)
        diff_mask(    2:la) = (diff_dv(2:la) > 0.0)    ! returns the logic index of the jumping grids
        diff_mask_1(1:la-1) = diff_mask(2:la)   ! one grid before the jumping grid    
                
        ! for example, Fella 2014 Fig 1
        ! if jump_ind = [6,10], then jump_ind_1 = [5,9]
        ! it means: deriv_Vhat(8) > deriv_Vhat(7)
        ! deriv_Vhat(14) > deriv_Vhat(13)
        ! deriv_Vhat(20) > deriv_Vhat(19)
        ! compute the boundaries        
                            
        if (count(diff_mask) > 0) then   ! non-concave region exists
            
            noncv_deriv_Vmax = maxval(deriv_Vhat(2:la),diff_mask(  2:la))   ! derivative of V ON     the jumping point, deriv_Vhat(10)
            noncv_deriv_Vmin = minval(deriv_Vhat(2:la),diff_mask_1(2:la))   ! derivative of V BEFORE the jumping point, deriv_Vhat(5)
            
            if (noncv_deriv_Vmax > deriv_Vhat(1)) then   ! Fella 2014 Fig 1, {a'2...a'11}
                noncv_ind_min = 1
            else
                ! choose minimal index from the range deriv_Vhat > noncv_deriv_Vmax deriv_Vhat(10) ---> [2]
                noncv_ind_min = minloc(deriv_Vhat, 1, deriv_Vhat > noncv_deriv_Vmax) + 1
            end if
                                
            if (noncv_deriv_Vmin < deriv_Vhat(la)) then
                noncv_ind_max = la
            else
                ! choose maximal index from the range deriv_Vhat < noncv_deriv_Vmin deriv_Vhat(5) ---> [11]
                noncv_ind_max = maxloc(deriv_Vhat, 1, deriv_Vhat < noncv_deriv_Vmin) - 1     
            end if
                                                                
        else   ! all region concave
                                
            noncv_ind_min = 0
            noncv_ind_max = 0
                                
        end if

        
        ! ===========================================================================================================================================
        ! global solution 
                            
        glob_ind(1:la) = 0  ! initialize global index                            
                            
        if (noncv_ind_min == 0 .and. noncv_ind_max == 0) then  ! no non-concave region exists
                                
            glob_ind(1:la) = 1
                                           
        else  ! non-concave region exists, noncv_ind_min >=1 or noncv_ind_max <= la
                               
            do ja = 1,la   
                                    
                if (ja < noncv_ind_min .or. ja > noncv_ind_max) then  ! lies outside the non-concave region; within concave region
                                        
                    glob_ind(ja) = 1
                                        
                else  ! lies inside the non-concave region, need further check to see if the global solution
                                        
                    coh_temp = coh(ja)  ! (1+(1-tau_a)*r)(a + bequest) + w*e*eps*l (sos) - taxw - taxss + lump 
                                        ! if global solution, then labor supply h and consumption c are also global
                                        ! if this cash-on-hand does not yield apgrid(ja) as a global solution
                                        ! then labor and consumption should also be modified
                    
                    c_temp = (coh_temp - apgrid(ja)) / (1.+tau_c)
                    Vhat_temp(1:la) = -1e10
                    ! the corresponding h_temp is indexed by h_vec(ja)
                    
                    if (c_temp > 0.) then
                        do jja = noncv_ind_max, noncv_ind_min, -1
                            ! search new a' from the highest apgrid(noncv_ind_max=11) to the lowest apgrid(noncv_ind_min=2) 
                            ! eventually stops at the lowerst a' such that c is positive
                            ! at this moment, c and u are maximized 
                            c_temp = (coh_temp - apgrid(jja)) / (1.+tau_c)
                            if (c_temp > 0) then
                                call sub_util(c_temp, h_vec(jja), u_temp, xx, xx)
                                Vhat_temp(jja) = u_temp + Vhat(jja)
                            else
                                Vhat_temp(jja) = -1e10
                            end if
                        end do
                    end if
                                      
                    ! for the newly computed Vhat_temp over the range of apgrid(noncv_ind_min) and apgrid(noncv_ind_max)
                    ! if the max Vhat_temp occurs within [ja-1, ja+1], then we consider the original index ja is a global solution
                    max_temp = maxval(Vhat_temp)                                        
                    ind_jja  = maxloc(Vhat_temp, dim = 1)
                                        
                    if (abs(ind_jja - ja) <= 1) then   ! allowing for one grid deviation, Youngsoo
                        glob_ind(ja) = 1
                    end if
                                                                                                                        
                end if
                
            end do  ! ja = 1:la                            
                                            
        end if   ! concave / non-concave regions judgements   
                        
        
        ! ===========================================================================================================================================
        ! find out the current agrid_pack, for interpolation later
                    
        do ja = 1,la  ! find the fist index of global solution
                                
            if (glob_ind(ja) == 1) then
                glob_ind_first = ja
                exit
            end if                                
                                
        end do
        

        ! compute packed apgrid, k and Vhat, 
        ! meaning on the packed set, all grids are global 
        
        apgrid_pack = 0.0
        a0_pack     = 0.0
        Vhat_pack   = 0.0
                            
        la_pack = count(glob_ind == 1)
        mask_glob_ind = glob_ind .eq. 1
                            
        if (noncv_ind_min == 0 .and. noncv_ind_max == 0) then  ! all concave region
            
            ! construct packed region directly
            apgrid_pack(1:la_pack) = pack(apgrid, mask_glob_ind)
            a0_pack(    1:la_pack) = pack(a0,     mask_glob_ind)
            Vhat_pack(  1:la_pack) = pack(Vhat,   mask_glob_ind)
            Vchat_pack( 1:la_pack) = pack(Vchat,  mask_glob_ind)
                                
        elseif (noncv_ind_min > 1) then  ! glob_ind_first == 1, apgrid(1) is the global solution, included in packed apgrid
                                
            ! construct packed region directly
            apgrid_pack(1:la_pack) = pack(apgrid, mask_glob_ind)
            a0_pack(    1:la_pack) = pack(a0,     mask_glob_ind)
            Vhat_pack(  1:la_pack) = pack(Vhat,   mask_glob_ind)
            Vchat_pack( 1:la_pack) = pack(Vchat,  mask_glob_ind)
                                
        elseif (noncv_ind_min == 1 .and. glob_ind_first == 1) then  ! non-concave region starts with 1st index, but apgrid(1) is a global solution
                            
            ! construct packed region directly
            apgrid_pack(1:la_pack) = pack(apgrid, mask_glob_ind)
            a0_pack(    1:la_pack) = pack(a0,     mask_glob_ind)
            Vhat_pack(  1:la_pack) = pack(Vhat,   mask_glob_ind)
            Vchat_pack( 1:la_pack) = pack(Vchat,  mask_glob_ind)
                                
        else  
                        
            ! noncv_ind_min == 1 .and. glob_ind_first > 1
            ! non-concave region starts from apgrid(1), which is not a global solution
            ! an additional ap_lower_bound should be inserted  
            ! construct temporary packed apgrids, will later add ap_lower_bound     
            apgrid_pack_temp(1:la_pack) = pack(apgrid, mask_glob_ind)
            a0_pack_temp(    1:la_pack) = pack(a0,     mask_glob_ind)
            Vhat_pack_temp(  1:la_pack) = pack(Vhat,   mask_glob_ind)
            Vchat_pack_temp( 1:la_pack) = pack(Vchat,  mask_glob_ind)
            

            ! a_lower_bound, or cutoff ap
            ! under this a_lower_bound, a' = apgrid(1), namely constrained
            ! Fella (2014), find a resource z such that
            ! u(z - apgrid(noncv_ind_min)) + Vhat(apgrid(noncv_ind_min)) = u(z - apgrid(glob_ind_first)) + Vhat(apgrid(glob_ind_first))
            ! e.g., noncv_ind_min = 1, glob_ind_first = 5
            ! the agrid corrosponding to the resource z is the a_lower_bound
            ! as resource increases, the RHS increases
            ! Thus, 
            ! step1: find a resource such that LHS > RHS, denote this resource resource1
            ! step2: gradually increases resource till LHS < RHS, denote this resource resource2
            ! step3: bisection to find a resource between resource1 and resource2 such that LHS exactly == RHS
            
            
            ! ---------------------------------------------------------------------------------------------------------------------------------------            
            ! guess a resource1, such that LHS > RHS, if not, gradually reduce resource 
            ! such a resource1 must yeild c > 0, if not, gradually increase resource1

            ! .......................................................................................................................................
            ! the way to pin down resource1, this method is too slow                                                                                !
            !                                                                                                                                       !
            ! resource1 = apgrid(1)                                                                                                                 !
            !                                                                                                                                       !
            ! c_noncv_ind_min = ( resource1 - apgrid(noncv_ind_min) ) / (1+tau_c)                                                                   !
            ! c_glob_first    = ( resource1 - apgrid(glob_ind_first)) / (1+tau_c)                                                                   !
            !                                                                                                                                       !
            ! do while (c_noncv_ind_min <= 0 .or. c_glob_first <= 0 )                                                                               !
            !                                                                                                                                       !
            !    resource1 = resource1 + 1e-8    ! + 1e-8 should be very small, such that c2_glob_first is epsilon-greater than 0                   !
            !                                     ! in order to guarantee u_glob_first is significantly small                                       !
            !                                                                                                                                       !
            !     c_noncv_ind_min = ( resource1 - apgrid(noncv_ind_min) ) / (1+tau_c)                                                               !
            !     c_glob_first    = ( resource1 - apgrid(glob_ind_first)) / (1+tau_c)                                                               !
            !                                                                                                                                       !
            ! end do                                                                                                                                !
            !........................................................................................................................................
            
            ! to get c_glob_first = 1e-16, in order to make sure u_glob_first is sufficiently small
            c_glob_first = 1e-16
            resource1 = (1+tau_c) * c_glob_first + apgrid(glob_ind_first)
            ! c_noncv_ind_min = ( resource1 - apgrid(noncv_ind_min) ) / (1+tau_c)
            !                 > ( resource1 - apgrid(glob_ind_first)) / (1+tau_c) = c_glob_first
                        
            ! now both c_noncv_ind_min > 0 and c_glob_first > 0 
            ! LHS = u(z - apgrid(noncv_ind_min)) + Vhat(apgrid(noncv_ind_min))
            c_noncv_ind_min = ( resource1 - apgrid(noncv_ind_min)) / (1+tau_c)
            
            ! uc_temp 
            call sub_util(c_noncv_ind_min, xx, xx, uc_temp, xx)
            
            ! corrosponding hrs 
            ! search a small region
            lh = 30
            hmin = max(0., min(h_vec(noncv_ind_min), h_vec(glob_ind_first) - 0.1)
            hmax = min(1. - 1e-8, max(h_vec(noncv_ind_min), h_vec(glob_ind_first) + 0.1)
            call sub_hrs(lh, hmin, hmax, c_noncv_ind_min, uc_temp, indh, h_temp)
                        
            call sub_util(c_noncv_ind_min, h_temp, u_noncv_ind_min, uc_noncv_ind_min, xx)
            V_noncv_ind_min  =  u_noncv_ind_min +  Vhat(noncv_ind_min)
            Vc_noncv_ind_min = uc_noncv_ind_min + Vchat(noncv_ind_min)
            
            
            ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! RHS = u(z - apgrid(glob_ind_first)) + Vhat(apgrid(glob_ind_first))
            c_glob_first = (resource1 - apgrid(glob_ind_first)) / (1+tau_c) 
            
            ! uc_temp 
            call sub_util(c_glob_first, xx, xx, uc_temp, xx)

            ! corrosponding hrs 
            ! search a small region
            lh = 30
            hmin = max(0., min(h_vec(noncv_ind_min), h_vec(glob_ind_first) - 0.1)
            hmax = min(1. - 1e-8, max(h_vec(noncv_ind_min), h_vec(glob_ind_first) + 0.1)
            call sub_hrs(lh, hmin, hmax, c_glob_first, uc_temp, indh, h_temp)
                        
            
            call sub_util(c_glob_first, u_glob_first, uc, un)
            V_glob_first = u_glob_first + Vhat(glob_ind_first)
            
            
            ! NOW: V_noncv_ind_min > V_glob_first
            if (V_noncv_ind_min < V_glob_first) then
                print*, "something's wrong, V_noncv_ind_min < V_glob_first"
            end if
            

            ! ---------------------------------------------------------------------------------------------------------------------------------------            
            ! guess a resource2, such that LHS < RHS, if not, gradually increase resource 
            ! assuming resource2 = resource1 + 1e-2, definitely such that c2_noncv_ind_min > 0 and c2_glob_first > 0 
            resource2 = resource1 + 1e-2 
            
            ! LHS = u(z - apgrid(noncv_ind_min)) + Vhat(apgrid(noncv_ind_min))
            c_noncv_ind_min = (resource2 - apgrid(noncv_ind_min)) / (1+tau_c)
            
            ! uc_temp 
            call sub_util(c_noncv_ind_min, xx, xx, uc_temp, xx)

            ! corrosponding hrs 
            ! search a small region
            lh = 30
            hmin = max(0., min(h_vec(noncv_ind_min), h_vec(glob_ind_first) - 0.1)
            hmax = min(1. - 1e-8, max(h_vec(noncv_ind_min), h_vec(glob_ind_first) + 0.1)
            call sub_hrs(lh, hmin, hmax, c_noncv_ind_min, uc_temp, indh, h_temp)
                                    
            call sub_util(c_noncv_ind_min, h_temp, u_noncv_ind_min, uc_noncv_ind_min, xx)
            V_noncv_ind_min  =  u_noncv_ind_min +  Vhat(noncv_ind_min)
            Vc_noncv_ind_min = uc_noncv_ind_min + Vchat(noncv_ind_min)
            
              
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! RHS = u(z - apgrid(glob_ind_first)) + Vhat(apgrid(glob_ind_first))
            c_glob_first = (resource2 - apgrid(glob_ind_first) ) / (1+tau_c)
            
            ! uc_temp 
            call sub_util(c_glob_first, xx, xx, uc_temp, xx)

            ! corrosponding hrs 
            ! search a small region
            lh = 30
            hmin = max(0., min(h_vec(noncv_ind_min), h_vec(glob_ind_first) - 0.1)
            hmax = min(1. - 1e-8, max(h_vec(noncv_ind_min), h_vec(glob_ind_first) + 0.1)
            call sub_hrs(lh, hmin, hmax, c_glob_first, uc_temp, indh, h_temp)
            
            call sub_util(c_glob_first, h_temp, u_glob_first, uc, xx)
            V_glob_first = u_glob_first + Vhat(glob_ind_first)
                        
            ! check if such resours2 makes V_noncv_ind_min < V_glob_first
            iter_v = 0
            do while (V_noncv_ind_min > V_glob_first)    ! till V_noncv_ind_min < V_glob_first
                
                iter_v = iter_v + 1
                if (iter_v > 1000) then
                    exit
                end if
                
                resource2 = resource2 + 1e-1
                
                ! LHS = u(z - apgrid(noncv_ind_min)) + Vhat(apgrid(noncv_ind_min))
                c_noncv_ind_min = (resource2 - apgrid(noncv_ind_min)) / (1+tau_c)
                
                ! uc_temp 
                call sub_util(c_noncv_ind_min, xx, xx, uc_temp, xx)

                ! corrosponding hrs 
                ! search a small region
                lh = 30
                hmin = max(0., min(h_vec(noncv_ind_min), h_vec(glob_ind_first) - 0.1)
                hmax = min(1. - 1e-8, max(h_vec(noncv_ind_min), h_vec(glob_ind_first) + 0.1)
                call sub_hrs(lh, hmin, hmax, c_noncv_ind_min, uc_temp, indh, h_temp)
                            
                call sub_util(c_noncv_ind_min, h_temp, u_noncv_ind_min, uc_noncv_ind_min, xx)
                V_noncv_ind_min  =  u_noncv_ind_min +  Vhat(noncv_ind_min)
                Vc_noncv_ind_min = uc_noncv_ind_min + Vchat(noncv_ind_min)
                          
        
                !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                ! RHS = u(z - apgrid(glob_ind_first)) + Vhat(apgrid(glob_ind_first))
                c_glob_first = (resource2 ) / (1+tau_c)
                
                ! uc_temp 
                call sub_util(c_glob_first, xx, xx, uc_temp, xx)

                ! corrosponding hrs 
                ! search a small region
                lh = 30
                hmin = max(0., min(h_vec(noncv_ind_min), h_vec(glob_ind_first) - 0.1)
                hmax = min(1. - 1e-8, max(h_vec(noncv_ind_min), h_vec(glob_ind_first) + 0.1)
                call sub_hrs(lh, hmin, hmax, c_glob_first, uc_temp, indh, h_temp)
                
                call sub_util(c_glob_first, h_temp, u_glob_first, uc_glob_first, xx)
                V_glob_first  = u_glob_first  +  Vhat(glob_ind_first) 
                Vc_glob_first = uc_glob_first + Vchat(glob_ind_first) 
                                
            end do
            
            ! NOW: V_noncv_ind_min < V_glob_first
            if (V_noncv_ind_min > V_glob_first) then
                print*, "something's wrong, V_noncv_ind_min > V_glob_first"
            end if
            
            ! ---------------------------------------------------------------------------------------------------------------------------------------            
            ! NOW resource1 such that V_noncv_ind_min > V_glob_first
            ! BUT resource2 such that V_noncv_ind_min < V_glob_first
            ! but both resource1 & resource2 are coarse resources
            ! bisection to refine resource
            erro_v = abs(V_noncv_ind_min - V_glob_first)
            iter_v = 0
            
            do while (erro_v > 1e-4)
                
                iter_v = iter_v + 1
                
                if (iter_v > 200) then
                    exit
                end if
                                
                 resource = 0.5 * (resource1 + resource2)
                
                ! LHS = u(z - apgrid(noncv_ind_min)) + Vhat(apgrid(noncv_ind_min))
                c_noncv_ind_min = (resource - apgrid(noncv_ind_min)) / (1+tau_c)
                
                ! uc_temp 
                call sub_util(c_noncv_ind_min, xx, xx, uc_temp, xx)

                ! corrosponding hrs 
                ! search a small region
                lh = 30
                hmin = max(0., min(h_vec(noncv_ind_min), h_vec(glob_ind_first) - 0.1)
                hmax = min(1. - 1e-8, max(h_vec(noncv_ind_min), h_vec(glob_ind_first) + 0.1)
                call sub_hrs(lh, hmin, hmax, c_noncv_ind_min, uc_temp, indh, h_temp)
                
                call sub_util(c_noncv_ind_min, h_temp, u_noncv_ind_min, uc_noncv_ind_min, xx)
                V_noncv_ind_min  =  u_noncv_ind_min +  Vhat(noncv_ind_min)
                Vc_noncv_ind_min = uc_noncv_ind_min + Vchat(noncv_ind_min)
        
                !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                ! RHS = u(z - apgrid(glob_ind_first)) + Vhat(apgrid(glob_ind_first))
                c_glob_first = (resource - apgrid(glob_ind_first) ) / (1+tau_c)
                
                ! uc_temp 
                call sub_util(c_glob_first, xx, xx, uc_temp, xx)

                ! corrosponding hrs 
                ! search a small region
                lh = 30
                hmin = max(0., min(h_vec(noncv_ind_min), h_vec(glob_ind_first) - 0.1)
                hmax = min(1. - 1e-8, max(h_vec(noncv_ind_min), h_vec(glob_ind_first) + 0.1)
                call sub_hrs(lh, hmin, hmax, c_glob_first, uc_temp, indh, h_temp)
                
                call sub_util(c_glob_first, h_temp, u_glob_first, uc_glob_first, un)
                V_glob_first  =  u_glob_first +  Vhat(glob_ind_first)
                Vc_glob_first = uc_glob_first + Vchat(glob_ind_first)
                
                ! update
                if(V_noncv_ind_min > V_glob_first) then
                    resource1 = resource
                else
                    resource2 = resource
                end if
                                
                erro_v = abs(V_noncv_ind_min - V_glob_first)
                
            end do
                        
            
            ! such resource makes next period capital  = 0, i.e., a' = apgrid(1)
            ! thus c1 = c1_noncv_ind_min, c2 = c1_noncv_ind_min, h = h_temp
            ! backing out the current cut-off a
            resource   = resource1
            if (age <= tr ) then
                yy_temp = w * eprof_age_ablt(age,n) * eps(is,age,n) * h_temp
                taxss_temp   = tau_ss * min(y_bar, yy_temp)
                call sub_tax_hsv(prog, yy_vec(ja) - 0.5* taxss_temp,taxw_temp, mrgtaxw_temp)                
            else
                yy_temp = replace * earn_ablt(n) 
                taxss_temp  = 0.
                call sub_tax_hsv(prog, yy_temp,taxw_temp, mrgtaxw_temp)                
            end if
                        
            a0_cutoff = ( resource - yy_temp + taxw_temp + taxss_temp -lump) / ( 1 + (1-tau_a)*r ) - bequest

                                  
            ! the yielded a0_cutoff is the k_min, such that k' = apgrid(1)
            ! adding a0_cutoff to the packed current k
            ! adding apgrid(1) to the packed k' (since k'(1) is not a globle solution, k'(1) was not in apgrid_pack_temp)
            apgrid_pack(1)             = apgrid(1)
            apgrid_pack(2:la_pack + 1) = apgrid_pack_temp(1:la_pack)
            a0_pack(1)                 = a0_cutoff
            a0_pack(2:la_pack + 1)     = a0_pack_temp(1:la_pack)                                                               
            Vhat_pack(1)               = u_noncv_ind_min + Vhat(noncv_ind_min)
            Vhat_pack(2:la_pack + 1)   = Vhat_pack_temp(1:la_pack)
            
            Vchat_pack(1)              = uc_noncv_ind_min + Vchat(noncv_ind_min)
            Vchat_pack(2:la_pack + 1)  = Vchat_pack_temp(1:la_pack)
                                
            la_pack = la_pack + 1  ! the length of global solution add 1
              
                                
        end if  ! packed a, a' & V      
        
                
        ! ===========================================================================================================================================
        ! interpolating x'_pack to the current x_pack
     
        do ia = 1,la  ! the current capital agrid(ia)  
            
            if (agrid(ia) <= a0_pack(1)) then 
                
                ! must binding, increase a0_pack(1) is not the minimum of a0_pack
                inda_vec(ia) = 1
                wgta_vec(ia) = 1.0
                a_pol_egm(ia) = 0.0
                                
            elseif (agrid(ia) >= a0_pack(la_pack)) then
                
                inda_vec(ia) = la_pack - 1
                wgta_vec(ia)  = 0.0
                a_pol_egm(ia) = apgrid(la_pack)

            else   ! a0_pack(1) < agrid(ia) < a0_pack(la)
                
                call sub_weight(la, a0_pack, agrid(ia), inda_vec(ia), wgta_vec(ia))
                a_pol_egm(ia) = wgta_vec(ia) * apgrid_pack(inda_vec(ia)) + (1.0 - wgta_vec(ia)) * apgrid_pack(inda_vec(ia)+1) 
            
            end if
            
            EVhat(ia) = wgta_vec(ia) * Vhat_pack( inda_vec(ia)) + (1.0 - wgta_vec(ia)) * Vhat_pack( inda_vec(ia)+1)
            EVchat(ia)= wgta_vec(ia) * Vchat_pack(inda_vec(ia)) + (1.0 - wgta_vec(ia)) * Vchat_pack(inda_vec(ia)+1)
            
        end do ! current capital ia

        ! retired 
        if (age > tr) then 
            do ia = 1,la  ! the current capital agrid(ia)
            
                ! B.C. (1+tau_c)*c + a_pol_egm(ia) = (1 + (1-tau_a)*r)*(agrid+bequest) + yy - yy_tax -yy_ss + lump 
                sos_pol_egm(ia) = replace * earn_ablt(n) 
                call sub_tax_hsv(prog, sos_pol_egm(ia), taxw_pol_egm(ia), xx)
                taxss_pol_egm(ia) = 0.
                
                h_pol_egm(ia) = 0.0
                
                ! resource
                resource = (1+(1-tau_a)*r) * (agrid(ia) + bequest) + sos_pol_egm(ia) - taxw_pol_egm(ia) - taxss_pol_egm(ia)  + lump
                
            end do
        
        else 
            ! workign hours worked and consumption are computed simulataneously
            do ia = 1,la
            
                if (ia == 1) then
                    lh = 100
                    hmin = 0.
                    hmax = 1-1e-8
                    
                    call sub_con_hrs(lh,hmin,hmax,a_pol_egm(ia),EVhat(ia), indh,h_temp)
                    
                else ! ia > 1 , hgrid is defined in a neighbourhood of h_pol_egm(ia-1)
                    
                    lh = 30
                    hmin = max(0.,     h_pol_egm(ia-1) - 0.1)
                    hmax = min(1-1e-8, h_pol_egm(ia-1) + 0.1)
                    
                    call sub_con_hrs(lh,hmin,hmax,a_pol_egm(ia),EVhat(ia), indh, h_temp)
                    
                    ! reducing lower bound
                    if (indh == 1 .and. hmin > 1e-8) then
                        lh = 30
                        hmin = 0.
                        hmax = h_pol_egm(ia-1) - 0.1
                        
                        call sub_con_hrs(lh,hmin,hmax,a_pol_egm(ia),EVhat(ia), indh, h_temp)
                        
                    end if
                    
                end if
                
                ! hours
                h_pol_egm(ia) = h_temp                
                yy_temp = w * eprof_age_ablt(age,n) * eps(is,age,n) * h_pol_egm(ia)
                taxss_pol_egm(ia)  = tau_ss * min(y_bar, yy_temp)
                call sub_tax_hsv(prog, yy_temp-0.5*taxss_pol_egm(ia),taxw_pol_egm(ia), xx)
                    
                ! resource
                resource = (1+(1-tau_a)*r) * (agrid(ia) + bequest) + yy_temp - taxw_pol_egm(ia) - taxss_pol_egm(ia) + lump 
                
            end do 
            
        end if ! current capital ia


        do ia = 1,la  ! the current capital agrid(ia)
            
            ! c policies
            c_pol_egm(ia) = ( resource - a_pol_egm(ia) ) / (1+tau_c) 
                
            if (c_pol_egm(ia) > 1e-4) then
                call sub_util(c_pol_egm(ia), h_pol_egm(ia), u_temp, uc_temp, xx)
                        
                V_pol_egm( ia) =  u_temp + EVhat( ia)
                Vc_pol_egm(ia) = uc_temp + EVchat(ia)
                Vn_pol_egm(ia) = V_pol_egm(ia) - Vc_pol_egm(ia)
                
            else
                
                V_pol_egm( ia) = -1e5
                Vc_pol_egm(ia) = -1e5
                Vn_pol_egm(ia) = -1e5
                
            end if                                                        
            
        end do ! current capital ia

        
    
end subroutine 
