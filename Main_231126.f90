!  Main_231126.f90 
!
!  FUNCTIONS:
!  Main_231126 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Main_231126
!
!  PURPOSE: calibration
!   
!  main programe loop: beta ---> data_KY; B ---> data_hrs; 
!  steady state loop : r; tau_c, lump ---> gov budget; replace ---> sos balance; bequest ---> balanced 
!  
!****************************************************************************

    program Main_231126
    
    use mod_gvar
    use mod_functions

    implicit none

    !====================================================================================
    ! exogenously determined parameters
    ! survival probability
    OPEN(3,file=fileplace_input//"/2y_survival_spline.txt",status='old',action='READ')
    do age = 1, maxage
        READ(3,*) phi_Age(age)
    end do
    
    ! percentage of each age
    cohsz_age(1) = 1
    do age = 2,maxage
        cohsz_age(age) = phi_age(age-1)*cohsz_age(age-1) / (1+g);
    end do
    poppct_age = cohsz_age / sum(cohsz_age)     
    
    !--------------------------------------------------------------------------
    ! initial distribution
    ! share of abilities
    pi_ablt(1) = 0.7
    pi_ablt(2) = 0.3
    
    !! call for approperiate shock
    !if (shockcase = 1 ) then
    !    ! tauchen
    !    ! output eps(ls,maxage,ln), pistar(ls,ln), pi(ls,ls,maxage,ln)
    !    
    !elseif (shockcase = 2 ) then
    !    ! conesa
    !    ! read eps(ls,maxage,ln), pi(ls,ls,maxage,ln)
    !    OPEN(3,file=fileplace_conesa//" .txt",status='old',action='READ')
    !    
    !    close(unit = 3)
    !    
    !    ! initial shock distribution
    !    age = 1
    !    pistar(1:ls,1) = (/1., 4., 5., 10., 20., 20., 20., 10., 5., 4., 1.  /) / 100.  
    !    pistar(1:ls,2) = pstar(1:ls,1) 
    !end if
    !
    !! initial distribution over assets, shock, ability
    !call sub_init_dist(la,ls,ln,prob_pol(la,ls,1,ln))
    prob_pol(1,1:ls,1,1) = (/1., 4., 5., 10., 20., 20., 20., 10., 5., 4., 1.  /) / 100. * pi_ablt(1)
    prob_pol(1,1:ls,1,2) = (/1., 4., 5., 10., 20., 20., 20., 10., 5., 4., 1.  /) / 100. * pi_ablt(2)
        
    
    !--------------------------------------------------------------------------
    ! discretize asset grid
    agrid(1) = 1.e-6
    agrid(la) = 200.0
    do ia = 2,la-1
        agrid(ia) = (sqrt(agrid(1)) + ia * (sqrt(agrid(la)) - sqrt(agrid(1))) / la)**2.0
    end do
    apgrid = agrid
    

    !--------------------------------------------------------------------------
    ! start calibration
    beta = 0.96
    B = 3
        
    

    end program Main_231126

