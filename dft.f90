module dft
    implicit none
    REAL(8), PARAMETER :: kB = 1.38064852d-23
    REAL(8), PARAMETER :: NA = 6.022140857d23
    REAL(8), PARAMETER :: pi = acos(-1.d0)
    REAL(8) :: Ta ! temperatura em kelvin
    REAL(8) :: T ! temperatura admensional
    REAL(8) :: P
    REAL(8) :: H
    REAL(8) :: sigma_ff 
    REAL(8) :: sigma_sf
    REAL(8) :: eps_ff
    REAL(8) :: eps_sf
    REAL(8) :: d
    REAL(8) :: rmin
    REAL(8) :: rc
    REAL(8) :: rho_s
    REAL(8) :: Delta
    REAL(8) :: dz
    REAL(8) :: tol
    
    contains 

    subroutine ler_dados()
        implicit none

        open(1, file='dados.inp', status='old')

        read(1,*) Ta 
        read(1,*) P 
        read(1,*) H 
        read(1,*) eps_ff 
        read(1,*) sigma_ff 
        read(1,*) eps_sf 
        read(1,*) sigma_sf 
        read(1,*) rho_s  
        read(1,*) Delta 
        read(1,*)
        read(1,*) dz
        read(1,*) tol
        
        close(1)

        P = (P*1.d5*(sigma_ff*1.d-10)**3)/(eps_ff*kB)
        T = Ta/eps_ff
        eps_sf = eps_sf/eps_ff
        eps_ff = eps_ff/eps_ff
        rho_s = rho_s*sigma_ff**3
        Delta = Delta/sigma_ff
        sigma_sf = sigma_sf/sigma_ff
        sigma_ff = sigma_ff/sigma_ff
        !d = sigma_ff
        d = sigma_ff*(0.3837d0*T+1.035d0)/(0.4249d0*T+1.d0)
        rmin = 2.d0**(1.d0/6.d0)*sigma_ff
        rc = 5.d0*sigma_ff
    
    end subroutine ler_dados

    function rtnewt(x0, xacc) 
        implicit none
        REAL(8), INTENT(IN) :: x0, xacc  
        REAL(8) :: rtnewt
        REAL(8) :: f, df, dx
    
        rtnewt = x0
        dx = 2.d0*xacc
        do while (abs(dx) > xacc)
            call funcd(rtnewt,f,df)
            dx = f/df
            rtnewt = rtnewt-dx
        end do
    
    end function rtnewt

    subroutine funcd(rho, fval, fderiv)
        implicit none
        REAL(8), INTENT(IN) :: rho
        REAL(8), INTENT(OUT) :: fval, fderiv
        REAL(8) :: cte, eta, Phi_att
        
        cte = 1.d0/6.d0
        eta = cte*(pi*rho*d**3)
        
        Phi_att = -dsqrt(2.d0)*(32.d0/9.d0)*pi*eps_ff*sigma_ff**3+(16.d0/3.d0)*pi*eps_ff*sigma_ff**3* &
        ((sigma_ff/rc)**3-(1.d0/3.d0)*(sigma_ff/rc)**9)

        fval = (rho*T)*((1.d0+eta+eta**2-eta**3)/(1.d0-eta)**3)+0.5d0*rho**2*Phi_att-P

        fderiv = 1.d0*Phi_att*rho + 3.d0*T*cte*d**3*pi*rho*(-cte**3*d**9*pi**3*rho**3 + cte**2*d**6*pi**2*rho**2 + &
        cte*d**3*pi*rho + 1.d0)/(-cte*d**3*pi*rho + 1.d0)**4 + T*rho*(-3.d0*cte**3*d**9*pi**3*rho**2 + & 
        2.d0*cte**2*d**6*pi**2*rho + cte*d**3*pi)/(-cte*d**3*pi*rho + 1.d0)**3 + T*(-cte**3*d**9*pi**3*rho**3 + & 
        cte**2*d**6*pi**2*rho**2 + cte*d**3*pi*rho + 1.d0)/(-cte*d**3*pi*rho + 1.d0)**3
    
    end subroutine funcd

    function chemical_potential(temp, dens) result(mu)
        implicit none
        REAL(8), INTENT(IN) :: temp, dens
        REAL(8) :: eta, mu, Phi_att 

        eta = (pi*dens*d**3)/6.d0 
        !mu = temp*(dlog(dens)+eta*(8.d0-9.d0*eta+3.d0*eta**2)/((1.0d0-eta)**3))
        mu = temp*(eta*(8.d0-9.d0*eta+3.d0*eta**2)/((1.d0-eta)**3))
        !print*, mu
        Phi_att = -dsqrt(2.d0)*(32.d0/9.d0)*pi*eps_ff*sigma_ff**3+(16.d0/3.d0)*pi*eps_ff*sigma_ff**3* &
        ((sigma_ff/rc)**3-(1.d0/3.d0)*(sigma_ff/rc)**9)
        !print*, dens*Phi_att 
        mu = mu+dens*Phi_att

    end function chemical_potential

    function trapz(x, y, a, b) result(r)
        implicit none
        REAL(8), INTENT(IN) :: x(:)
        REAL(8), INTENT(IN) :: y(:)
        INTEGER, INTENT(IN) :: a
        INTEGER, INTENT(IN) :: b
        REAL(8) :: r

        associate(n => size(x))
            r = 0.5d0*sum((y(a+1:b-0)+y(a+0:b-1))*(x(a+1:b-0)-x(a+0:b-1))) 
        end associate 

    end function trapz

    subroutine WD(kz, z, rho_in, ni)
        implicit none
        INTEGER, INTENT(IN) :: kz
        REAL(8), DIMENSION(:), INTENT(IN) :: z, rho_in
        REAL(8), DIMENSION(0:5), INTENT(OUT) :: ni
        REAL(8) :: zin, zsu 
        INTEGER :: kin, ksu 
    
        zin = z(kz)-0.5d0*d
        zsu = z(kz)+0.5d0*d
        kin = max(idnint(zin/dz), 1)
        ksu = min(idnint(zsu/dz), size(z))
    
        if ( kin < ksu ) then
            ni(2) = pi*d*trapz(z, rho_in, kin, ksu)
            ni(1) = ni(2)/(2.d0*pi*d)
            ni(0) = 2.d0*ni(1)/d
            ni(3) = pi*trapz(z, rho_in*(0.25d0*d**2-(z-z(kz))**2), kin, ksu)  
            ni(5) = -2.d0*pi*trapz(z, rho_in*(z-z(kz)), kin, ksu)   
            ni(4) = ni(5)/(2.d0*pi*d)    
        end if
    
    end subroutine WD

    function DA_HS(kz, z, rho_in) result(DA)
        implicit none
        INTEGER, INTENT(IN) :: kz
        REAL(8), DIMENSION(:), INTENT(IN) :: z, rho_in
        REAL(8) :: DA
        REAL(8) :: zin, zsu, cte
        REAL(8), DIMENSION(size(z)) :: dPhi_n0, dPhi_n1, dPhi_n2, dPhi_n3, dPhi_n4, dPhi_n5
        INTEGER :: i, n, kin, ksu
        REAL(8), DIMENSION(0:5) :: Ni

        n = size(z)

        do i = 1, n
            call WD(i, z, rho_in, Ni)
            cte = 1.d0-Ni(3) 
            dPhi_n0(i) = -dlog(cte)
            dPhi_n1(i) = Ni(2)/(cte)
            dPhi_n2(i) = Ni(1)/cte+3.d0*(Ni(2)**2-Ni(5)**2)*(Ni(3)+cte**2*dlog(cte))/(36.d0*pi*Ni(3)**2*cte**2)  
            dPhi_n3(i) = Ni(0)/cte+(Ni(1)*Ni(2)-Ni(4)*Ni(5))/(cte**2)-(Ni(2)**3-3.d0*Ni(2)*Ni(5)**2) &
            *(Ni(3)*(Ni(3)**2-5.d0*Ni(3)+2.d0)+2.d0*cte**3*dlog(cte))/(36.d0*pi*Ni(3)**3*cte**3)
            dPhi_n4(i) = -(Ni(5)/cte)
            dPhi_n5(i) = -(Ni(4)/cte)-Ni(2)*Ni(5)*(Ni(3)+cte**2*dlog(cte))/(6.d0*pi*Ni(3)**2*cte**2)
        end do
        
        zin = z(kz)-0.5d0*d
        zsu = z(kz)+0.5d0*d
        kin = max(idnint(zin/dz), 1)
        ksu = min(idnint(zsu/dz), size(z))

        if ( kin < ksu ) then
            DA = trapz(z, (1.d0/d)*dPhi_n0, kin, ksu)
            DA = DA + trapz(z, 0.5d0*dPhi_n1, kin, ksu)
            DA = DA + trapz(z, pi*d*dPhi_n2, kin, ksu)
            DA = DA + trapz(z, pi*dPhi_n3*(0.25d0*d**2-(z-z(kz))**2), kin, ksu)
            DA = DA + trapz(z, (1.d0/d)*dPhi_n4*(z-z(kz)), kin, ksu)
            DA = DA + trapz(z, 2.d0*pi*dPhi_n5*(z-z(kz)), kin, ksu)
            DA = DA*T
        end if

    end function DA_HS

    function DA_DISP(kz, z, rho_in) result(DA)
        implicit none
        INTEGER, INTENT(IN) :: kz
        REAL(8), DIMENSION(:), INTENT(IN) :: z, rho_in
        REAL(8) :: DA
        REAL(8), DIMENSION(size(z)) :: DAv
        REAL(8), DIMENSION(:), ALLOCATABLE :: DAf
        INTEGER :: i, j, n, n1, n2

        n = size(z)
        n1 = 0
        n2 = 0
        DAv = 0.d0
        do i = 1, n
            if ( abs(z(i)-z(kz)) > rmin .and. abs(z(i)-z(kz)) < rc ) then
                DAv(i) = (0.4d0*eps_ff*sigma_ff**12)*(1.d0/abs(z(i)-z(kz))**10-1.d0/rc**10)+ &
                (eps_ff*sigma_ff**6)*(1.d0/rc**4-1.d0/abs(z(i)-z(kz))**4)
                DAv(i) = rho_in(i)*DAv(i)
                n1 = n1+1
            else if ( abs(z(i)-z(kz)) < rmin ) then
                DAv(i) = 0.5d0*eps_ff*(abs(z(i)-z(kz))**2-rmin**2)+(0.4d0*eps_ff*sigma_ff**12)*(1.d0/rmin**10-1.d0/rc**10)+ &
                (eps_ff*sigma_ff**6)*(1.d0/rc**4-1.d0/rmin**4)
                DAv(i) = rho_in(i)*DAv(i)
                n2 = n2+1
            end if
        end do

        ALLOCATE(DAf(n1+n2))

        j = 1
        do i = 1, n
            if ( DAv(i) /= 0.d0 ) then
                DAf(j) = DAv(i)
                j = j+1
            end if
        end do

        DA = 2.d0*pi*trapz(z, DAf, 1, n1+n2)

        DEALLOCATE(DAf)
        
    end function DA_DISP

    function external_potential(z) result(V)
        implicit none
        REAL(8), INTENT(IN) :: z
        REAL(8) :: V
        REAL(8) :: parametro

        parametro = 2.d0*pi*rho_s*eps_sf*sigma_sf**2*Delta
        V = parametro*(0.4d0*(sigma_sf/z)**10-(sigma_sf/z)**4-sigma_sf**4/(3.d0*Delta*(z+0.61d0*Delta)**3))
        V = V + parametro*(0.4d0*(sigma_sf/(H-z))**10-(sigma_sf/(H-z))**4-sigma_sf**4/(3.d0*Delta*((H-z)+0.61d0*Delta)**3))
        if ( V > 100.d0 .or. isnan(V) ) then
            V = 100.d0
        end if
        !V = 10.d0*((1.d0/15.d0)*(sigma_sf/z)**9 - 0.5d0*(sigma_sf/z)**3)
        !V = V + 10.d0*((1.d0/15.d0)*(sigma_sf/(H-z))**9 - 0.5d0*(sigma_sf/(H-z))**3)

    end function external_potential

    subroutine picard(z, rho_in, rho_b, mu, Vext, rho_out)
        implicit none
        REAL(8), INTENT(IN) :: rho_b, mu
        REAL(8) :: alpha, soma_erro
        INTEGER :: i, it, Np
        REAL(8), DIMENSION(:), INTENT(IN) :: z, Vext
        REAL(8), DIMENSION(:), INTENT(INOUT) :: rho_in
        REAL(8), DIMENSION(:), INTENT(OUT) :: rho_out
        REAL(8), DIMENSION(size(z)) :: erro

        Np = size(z)
        soma_erro = 2.d0*tol
        it = 0
        alpha = 1.d-3

        do while (soma_erro > tol)
            do i = 1, Np
                rho_out(i) = rho_b*dexp((1.d0/T)*(mu-DA_HS(i, z, rho_in)-DA_DISP(i, z, rho_in)-Vext(i)))
                rho_out(i) = rho_in(i)*(1.d0-alpha)+rho_out(i)*alpha 
                erro(i) = abs(rho_out(i)-rho_in(i)) 
                rho_in(i) = rho_out(i)
            end do
            soma_erro = sum(erro)
            it = it+1
            alpha = 1.d-1
            print*, it, soma_erro 
        end do
        
    end subroutine picard

end module dft
!
! \\ Programa principal \\
!
program main_dft
    use dft
    implicit none 
    REAL(8) :: rho_b, mu
    REAL(8) :: start, finish
    INTEGER :: i, Np
    REAL(8), DIMENSION(:), ALLOCATABLE :: z, Vext  
    REAL(8), DIMENSION(:), ALLOCATABLE :: rho_in, rho_out
    
    call cpu_time(start)

    call ler_dados

    Np = idnint(H/dz)-1
    
    ALLOCATE(z(Np))
    ALLOCATE(Vext(Np))
    ALLOCATE(rho_in(Np))
    ALLOCATE(rho_out(Np))

    rho_b = rtnewt(P/T, 1d-8)
    !print*, rho_b/(3.575d-10)**3*(1.d0/NA)

    !print*, P

    mu = chemical_potential(T, rho_b)
    
    do i = 1, Np
        z(i) = dble(i)*dz
        Vext(i) = external_potential(z(i))
    end do

    rho_in = rho_b
    
    call picard(z, rho_in, rho_b, mu, Vext, rho_out)

    open(2, file='saida.dat', status='replace')
    do i = 1, Np
        write(2,*) z(i), rho_out(i)
    end do
    
    close(2)

    call cpu_time(finish)
    print*, 'Time (s):', finish-start

end program main_dft
