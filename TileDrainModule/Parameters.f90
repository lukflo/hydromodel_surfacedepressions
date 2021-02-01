module global 

    implicit none 
    ! Variables
    real :: g = 9.81, pi = 3.14159265359
    double precision :: delt, delt1, k, olddelt
    integer :: segs, stressperiods
    
    double precision, dimension(:), allocatable :: delx, theta, L, D, mann, So, tol_y, Loca, tol_trib, tol_trib1
    double precision, dimension(:), allocatable :: tol_y1, tol_Q, tol_Q1
    integer, dimension(:), allocatable :: chantype, RecSeg, DSBCType, USBCType, N
    
    integer :: i, days! time index, space index, number of reaches, number of stressperiods
    integer :: id, flag, r, jj, converge, unitql, last_trib, jjj, iii
    integer :: iter_count, solv_count, loop_count    

    double precision :: alpha, slot_width, drain_coeff, exp_coeff, depth_to_drain, ang2
    
    double precision, dimension(:,:,:), allocatable :: Q, y, qtrib
    double precision, dimension(:,:), allocatable :: wt_day, Qlast, ylast
    double precision, dimension(:), allocatable :: ql, wt
    
    double precision, dimension(:), allocatable :: reach_off, resid, x, guessQ, guessy
    double precision, dimension(:,:), allocatable :: Je, J, IC_Q, IC_y, Q_up, Q_upbase,y_down
    double precision, dimension(:), allocatable :: oldQ, oldy
    
    
    real :: E, El, Er, decay
    real, dimension(:,:,:), allocatable :: phi, ctrib
    double precision, dimension(:,:), allocatable :: GW_Conc, GW_Load, Conc_BC
    integer :: transport_flag, restart_flag
    
    
    logical :: exist
    
end module global 