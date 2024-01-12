module mod_gvar
    
    ! folders
    CHARACTER(*), PARAMETER :: fileplace_main   = 'E:\AlienWork\ConesaLiWang_Shock\Main_231126', &
                            &  fileplace_conesa = 'E:\AlienWork\ConesaLiWang_Shock\Main_231126_conesa', &
                            &  fileplace_input  = 'E:\AlienWork\ConesaLiWang_Shock\Main_231126\input'
    
    ! integers
    integer, PARAMETER :: shockcase = 1, & ! 1:tauchen; 2:conesa
                        ! la: size of asset grid: agrid; ls = size of shock grid: eps; ln: education level; maxage
                        & la = 200, ls = 11, ln = 2, maxage = 60, tr = 45
    
    double precision, PARAMETER :: A = 1, alpha = 0.36, delta = 0.069, sigma = 2, xi = -0.5, & ! A: tfp; delta: annual depreciation rate; sigma: CRRA; xi: labor elasticity
                                ! population growth rate g
                                &   g = 0.012, &
                                ! government
                                ! progs: 0 flat tax, 1 progressive tax; tau_w: labor income tax; lambda_l HSV tax parameter
                                ! tau_k: capital tax; 
                                ! mu_lump: share of lump-sum transfer in gov revenue
                                ! mu_govc: share of gov consumption in gov revenue
                                &   prog = 0, tau_w = 0.27, lambda_l = 0.1, tau_a = 0.39, mu_lump = 0.04, mu_govc = 0.2, & 
                                ! calibration targets
                                ! Hrs: hours worked; ky: capital to output ratio; govy: government to output ratios; debty: debt to gdp ratio
                                &   data_hrs = 0.3, data_KY = 3, data_govy = 0.2, data_debty = 0.6, &
                                ! precision threshold
                                &   threshold_ss = 1e-4, &
                                ! updating speed
                                &   update = 0.1
    
    ! integers
    integer :: ia, is, n, age ! ia: index of agrid; is: index of eps; n: index of edu;
    
    ! real
    double precision :: agrid(la), apgrid(la), eprof(maxage, ln), eps(ls,maxage,ln), pi(ls,ls,maxage,ln), pistar(ls,ln), & 
                    ! agrid: asset grid; eprof: life cycle profile; eps: shock grid; pi: transition matrix; pistar: shock distribution at age 1
                    ! life-cycle
                    ! phi_age: age-dependent surival prob; cohsz_age: cohort size; poppct_age: population percentage by age
                    &   phi_age(maxage), cohsz_age(maxage), poppct_age(maxage), eprof_age_ablt(maxage,ln), &
                    ! preference
                    ! beta, B: disutility of work
                    &   beta, B, & 
                    ! prices: interest rate, wage rate
                    &   r, w, & 
                    ! ernings
                    &   y_bar, &
                    ! social security
                    ! tau_c: consumption tax; tau_ss: social security tax; replace: replacement ratio; ss: social security benefits
                    ! lump: lump-sum transfer; 
                    &   tau_c, tau_ss, replace, sos, lump, bequest, &
        
                    ! individual policy function
                    !  a_pol: policy funciton for assets; h_pol: hours worked; pol_l: effective labor supply; c_pol: consumption
                    !  v_pol: value; vc_pol: consumption part of value; vl_pol: labor part of value; 
                    !  taxc_pol: consumption tax payments; taxw_pol: labor tax payments; taxa_pol: capital tax payments; taxss_pol(la,ls,maxage,ln): 
                    !  social security tax payments; sos_pol: social security benefits                    
                    &   a_pol(la,ls,maxage,ln), h_pol(la,ls,maxage,ln), el_pol(la,ls,maxage,ln), c_pol(la,ls,maxage,ln), &
                    &   v_pol(la,ls,maxage,ln), vc_pol(la,ls,maxage,ln), vl_pol(la,ls,maxage,ln), &
                    &   taxc_pol(la,ls,maxage,ln), taxw_pol(la,ls,maxage,ln), taxa_pol(la,ls,maxage,ln), taxss_pol(la,ls,maxage,ln), sos_pol(la,ls,maxage,ln), &
                    ! distribution
                    ! sum(prob_pol(:,:,age,n)) = 1; sum(prob_unit_pol) = 1
                    &   prob_pol(la,ls,maxage,ln), prob_unit_pol(la,ls,maxage,ln), &
        
                    ! (maxage,ln)
                    &   bequest_age_ablt(maxage,ln), Kap_age_ablt(maxage,ln), hrs_age_ablt(maxage,ln), eLab_age_ablt(maxage,ln), earn_age_ablt(maxage,ln), &
                    &   Con_age_ablt(maxage,ln), Wel_age_ablt(maxage,ln), Welc_age_ablt(maxage,ln), Well_age_ablt(maxage,ln), &
                    &   taxa_age_ablt(maxage,ln), taxw_age_ablt(maxage,ln), taxc_age_ablt(maxage,ln), taxss__age_ablt(maxage,ln), sos_age_ablt(maxage,ln), &
    
                    ! sum up by abiltiy
                    ! pi_ablt: ability percentage
                    &   pi_ablt(ln), bequest_ablt(ln), Kap_ablt(ln), hrs_ablt(ln), ELab_ablt(ln), earn_ablt(ln), &
                    &   Con_ablt(ln), Wel_ablt(ln), Welc_ablt(ln), Well_ablt(ln), &
                    &   taxa_ablt(ln), taxw_ablt(ln), taxc_ablt(ln), taxss_ablt(ln), sos_ablt(ln), &
        
                    ! by age
                    &   kap_age(maxage), &
    
                    ! aggregates
                    ! kap_mod: capital; elab_mod: effective labor; hrs_mod: hours worked; earn_mod: earning; con_mod: consumption
                    ! V_mod: welfare; Vc_mod: consumption part of welfare; Vl_mod: labor part of welfare
                    &   gdp_mod, kap_mod, hrs_mod, hrs_wrk_mod, elab_mod, elab_wrk_mod, earn_wrk_mod, con_mod, V_mod, Vc_mod, Vl_mod, &
                    ! total population
                    &   pop_mod, &
                    ! tax
                    &   taxa_mod, taxw_mod, taxc_mod, taxss_mod, sos_mod, &
                    ! government
                    ! gov_mod: tax revenue; debt_mod: debt; lump: lump-sum transfer
                    &   gov_mod, debt_mod, &
                    ! ratios
                    &   ky_mod, govy_mod, govc_mod, debty_mod
        
                    
    
    
    
end module 