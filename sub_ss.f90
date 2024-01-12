!==============================================================================
!
! this subroutine conputes the steady state, given calibrated parameters
! iterating on: r; tau_c, lump ---> gov budget; replace ---> sos balance; bequest ---> balanced 
! 
!==============================================================================
    
    subroutine sub_ss
    
    use mod_gvar   
    use mod_functions
    
    implicit none  
    
    integer :: iter_r, tt, js, ja, ind_a
    
    double precision :: r_s, replace_s, bequest_s, lump_s, tau_cs, & 
                    ! prices etc to be iterated on
                    ! hrs worked grid
                    &   hgrid(ls), &
                    ! error term
                    &   erro_ss, &
                    ! utility
                    &   u, uc, ul, V_next(la), Vc_next(la), Vl_next(la), &
                    !
                    &   wgt_a
    
    ! initial guess of aggregates and prices
    GDP_mod         = 1.
    gov_mod         =  data_govy * GDP_mod
    earn_ablt(1:ln) = (1.-alpha) * GDP_mod
    
    r       = 0.04
    replace = 0.6
    bequest = 0.0
    lump    = mu_lump * gov_mod
    tau_c   = 0.05
    
    r_s       = r
    replace_s = replace 
    bequest_s = bequest
    lump_s    = lump 
    tau_cs    = tau_c
    
    erro_ss = 1
    
    do while (erro_ss > threshold_ss) 
        
        ! update prices
        r = (1.-update) * r + update * r_s         
        
        ! wage
        w = A*(1 - alpha)*((delta + r)/(A*alpha))**(alpha/(alpha - 1))
        
        ! new 
    
        ! initialization
        if ( iter_r == 1) then
            
            do n = 1,ln
                do age = 1,maxage
                    do is = 1,ls
                        a_pol(1:la,is,age,n) = agrid(1:la)
                    end do
                end do
            end do
            
        end if
        
        c_pol  = 0.0
        h_pol  = 0.0
        el_pol = 0.0
            
        V_pol   = -1e10
        Vc_pol  = -1e10
            
        taxc_pol  = 0.0
        taxw_pol  = 0.0
        taxa_pol  = 0.0
        taxss_pol = 0.0
            
        sos_pol = 0.0
            
        prob_pol(1:la,1:ls,2:maxage,1:ln) = 0.0
        
        !========================================================================================================================
        ! STEP1: policies on agrid, given all prices       
        ! last period
        age = maxage
        
        ! labor policy
        h_pol(1:la,1:ls,age,1:ln) = 0.0
            
        do n = 1,ln ! ability index, deterministic
                                                                                                    
            ! labor tax on social security benefits, depending on ability
            ia = 1
            is = 1
            sos_pol(1,1,age,n) = replace * earn_ablt(n)
            call sub_tax_hsv(sos_pol(1,1,age,n), taxw_pol(1,1,age,n))
            
            ! socail security benifits and tax
            sos_pol( 1:la,1:ls,age,n) = sos_pol( 1,1,age,n)
            taxw_pol(1:la,1:ls,age,n) = taxw_pol(1,1,age,n)
                    
            ! tax payments            
            is = 1
            do ia = 1,la

                ! capital
                a_pol(ia,is,age,n) = 0.0
                
                ! tax 
                taxc_pol( ia,is,age,n) = tau_c * c_pol(ia,is,age,n)
                taxa_pol( ia,is,age,n) = tau_a * a_pol(ia,is,age,n)
                taxss_pol(ia,is,age,n) = 0.
                    
                ! c
                c_pol(ia,is,age,n) = ((1 + (1-tau_a)*r) * agrid(ia) + sos_pol(ia,is,age,n) - taxw_pol(ia,is,age,n) - a_pol(ia,is,age,n)) / (1. + tau_c)
                    
                ! utility
                call sub_util(c_pol(ia,is,age,n), h_pol(ia,is,age,n), u, uc, ul)
                
                v_pol( ia,is,age,n) = u
                vc_pol(ia,is,age,n) = uc
                vl_pol(ia,is,age,n) = ul
                
                                
            end do   
            
            do is = 2,ls
                a_pol(    1:la,is,age,1:ln) = a_pol(    1:la,1,age,1:ln)
                c_pol(    1:la,is,age,1:ln) = c_pol(    1:la,1,age,1:ln)
                v_pol(    1:la,is,age,1:ln) = v_pol(    1:la,1,age,1:ln)
                vc_pol(   1:la,is,age,1:ln) = vc_pol(   1:la,1,age,1:ln)
                vl_pol(   1:la,is,age,1:ln) = vl_pol(   1:la,1,age,1:ln)
                taxc_pol( 1:la,is,age,1:ln) = taxc_pol( 1:la,1,age,1:ln)
                taxa_pol( 1:la,is,age,1:ln) = taxa_pol( 1:la,1,age,1:ln)
                taxss_pol(1:la,is,age,1:ln) = taxss_pol(1:la,1,age,1:ln)
            end do
            
        end do  !n
        
        
        ! -----------------------------------------------------------------------------------------------------------------------
        ! working and retired 
        do n = 1,ln
            
            do tt = 1,maxage-1
                age = maxage - tt ! from maxage-1:1
                
                do is = 1,ls
                    
                    do ja = 1,la               
                        V_next( ja) = beta * phi_age(age) * dot_product(pi(is,:,age,n), V_pol( ja,1:ls,age+1,n))          
                        Vc_next(ja) = beta * phi_age(age) * dot_product(pi(is,:,age,n), Vc_pol(ja,1:ls,age+1,n)) 
                        !Vl_next(ja) = beta * phi_age(age) * dot_product(pi(is,:,age,n), Vl_pol(ja,1:ls,age+1,n)) 
                    end do  
                
                    ! call EGM loop                            
                    call sub_EGM(   V_next(1:la),           Vc_next(1:la)                                  , &
                            &        a_pol(1:la,is,age,n),    c_pol(1:la,is,age,n),    h_pol(1:la,is,age,n), &
                            &        V_pol(1:la,is,age,n),   vc_pol(1:la,is,age,n)                         , &
                            &     taxc_pol(1:la,is,age,n), taxw_pol(1:la,is,age,n), taxa_pol(1:la,is,age,n), &
                            &    taxss_pol(1:la,is,age,n),  sos_pol(1:la,is,age,n)) 
                
                end do
                
                if (age <= tr) then
                    do is = 1,ls
                        do ia = 1,la
                            el_pol(ia,is,age,n) = eprof_age_ablt(age,n) * eps(is,age,n) * h_pol(ia,is,age,n)
                        end do
                    end do
                end if
                    
            end do ! age
        end do     ! n
                    
        
        !========================================================================================================================
        ! STEP2: distribution
        ! big transition matrix: prob((a,s) ----> (a',s')) or prob((a',s')|(a,s))
        ! without creating the big transition matrix
        ! prob(a',s') = sum_{a,s}[prob(a,s)*pi(a,a')*pi(s,s')]
        ! where pi(s,s') is the transition matrix for         
        ! initialize already 
        
        do n = 1,ln            
            do age = 1,maxage-1

                do is = 1,ls
                    do ia = 1,la      ! current asset and labor
                                            
                        ! find weight of a_pol(ia,is,age,n) on agrid
                        ! return wgt_a, the weight to ia; thus 1-wgt_a is the weight to ia+1
                        call sub_weight(la, apgrid, a_pol(ia,is,age,n), ind_a, wgt_a)
                             
                        ! goes from ia  to ind_a with probability of wgt_a
                        ! goes from is  to js    with probability of pi(is,js,age,n)
                        ! goes from age to age+1 with probabiltiy of 1.0
                        ! stay at   n            with probabiltiy of 1.0
                        do js = 1,ls
                            prob_pol(ind_a,  js,age+1,n) = prob_pol(ind_a,  js,age+1,n) + prob_pol(ia,is,age,n) * pi(is,js,age,n) *      wgt_a 
                            prob_pol(ind_a+1,js,age+1,n) = prob_pol(ind_a+1,js,age+1,n) + prob_pol(ia,is,age,n) * pi(is,js,age,n) * (1.0-wgt_a)
                        end do
                        
                    end do  ! ia
                end do      !is                 
                
                if ( abs(sum(prob_pol(:,:,age+1,n)) - 1.0) > 1e-6 ) then
                    print*, '       ----> something wrong with sum(prob_pol(:,:,age,n)):'
                    print*, n, age, sum(prob_pol(:,:,age+1,n))   
                end if  
                
                ! normalization by age and ability
                prob_pol(:,:,age+1,n) = prob_pol(:,:,age+1,n) / sum(prob_pol(:,:,age+1,n))
                
            end do ! age
            
        end do  ! skill
        
        
        ! whole population normalization
        do n = 1,ln
            do age = 1,maxage
                
                prob_unit_pol(:,:,age,n)   = prob_pol(:,:,age,n) * poppct_age(age) * pi_ablt(n) 
                                   
            end do ! age till 15
        end do

        if ( abs(sum(prob_unit_pol) - 1.0) > 1e-4 ) then
            print*, '       ----> something wrong with sum(prob_unit_pol):', sum(prob_unit_pol)
        end if
        
        prob_unit_pol = prob_unit_pol / sum(prob_unit_pol)
               
        
        !========================================================================================================================
        ! STEP3: aggregation

        ! (maxage,ln) :: 
        ! bequest_age_ablt, Kap_age_ablt, hrs_age_ablt, eLab_age_ablt, earn_age_ablt
        ! Con_age_ablt, Wel_age_ablt, Welc_age_ablt, Well_age_ablt, 
        ! taxa_age_ablt, taxw_age_ablt, taxc_age_ablt, taxss__age_ablt, sos_age_ablt
        
         ! (ln)::
         ! bequest_ablt, Kap_ablt, hrs_ablt, ELab_ablt, earn_ablt
         ! Con_ablt, Wel_ablt, Welc_ablt, Well_ablt
         ! taxa_ablt, taxw_ablt, taxc_ablt, taxss_ablt, sos_ablt
        
         ! (1) ::
         ! Kap, hrs_wrk, ELab, elab_wrk, earn_wrk
         ! Con, V, Vc, Vl
         ! taxa, taxw, taxc, taxss, sos
         ! gdp, KY, gov, govc, govy, debt, debty
                   
        ! aggregage 
        call sub_aggregation(la,ls,maxage,ln,    a_pol, prob_unit_pol,  kap_mod)
        call sub_aggregation(la,ls,maxage,ln,    c_pol, prob_unit_pol,  con_mod)
        call sub_aggregation(la,ls,maxage,ln,   el_pol, prob_unit_pol, elab_mod)
        call sub_aggregation(la,ls,maxage,ln, taxc_pol, prob_unit_pol, taxc_mod)
        call sub_aggregation(la,ls,maxage,ln, taxa_pol, prob_unit_pol, taxa_mod)
        call sub_aggregation(la,ls,maxage,ln, taxw_pol, prob_unit_pol, taxw_mod)
        
        call sub_aggregation(la,ls,tr,ln,    h_pol, prob_unit_pol(1:la,1:ls,1:tr,1:ln),  hrs_wrk_mod)
        call sub_aggregation(la,ls,tr,ln,   el_pol, prob_unit_pol(1:la,1:ls,1:tr,1:ln), elab_wrk_mod)
        call sub_aggregation(la,ls,tr,ln,taxss_pol, prob_unit_pol(1:la,1:ls,1:tr,1:ln),taxss_mod)
                
        call sub_aggregation(la,ls,maxage-tr,ln,sos_pol, prob_unit_pol(1:la,1:ls,tr+1:maxage,1:ln),sos_mod)
        
        call sub_aggregation(la,ls,maxage,ln,  V_pol, prob_unit_pol,  V_mod)
        call sub_aggregation(la,ls,maxage,ln, Vc_pol, prob_unit_pol, Vc_mod)
        call sub_aggregation(la,ls,maxage,ln, Vl_pol, prob_unit_pol, Vl_mod)
        
        earn_wrk_mod = w * elab_wrk_mod
        gdp_mod = A * Kap_mod**alpha * ELab_mod**(1-alpha)
        ky_mod = kap_mod / gdp_mod
        gov_mod = taxc_mod + taxa_mod + taxw_mod
        govy_mod = gov_mod / gdp_mod
        govc_mod = mu_govc * gov_mod
        debt_mod = (gov_mod - lump - govc_mod) / r
        debty_mod = debt_mod / gdp_mod
        
        !========================================================================================================================
        ! STEP4: updating prices
        ! interest rate
        r_s = A * alpha * Kap_mod**(alpha-1) * ELab_mod**(1-alpha) - delta
                
        ! pension
        replace_s = replace * (taxss_mod / sos_mod)**(0.1)
        
        ! bequest
        pop_mod = dot_product(poppct_age, phi_age) 
        bequest_s = bequest / pop_mod / (1.0+g) - poppct_age(1) * kap_age(1)
        
        ! lump-sum transfer
        lump_s = mu_lump * gov_mod
        
        ! tau_c
        ! gov budget: gov = r*debt + lump + govc
        tau_cs = (r*debt_mod + lump + govc_mod - taxw_mod - taxa_mod) / con_mod
        
        ! error
        erro_ss = max(abs(r-r_s), abs(replace_s - replace), abs(bequest_s - bequest), abs(lump_s-lump), abs(tau_cs - tau_c))
    
    end do
    
end subroutine