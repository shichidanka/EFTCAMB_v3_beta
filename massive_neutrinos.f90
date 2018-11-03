    module MassiveNu
    use constants, only: dl, const_pi
    implicit none
    private

    real(dl), parameter  :: fermi_dirac_const  = 7._dl/120*const_pi**4 ! 5.68219698_dl
    !fermi_dirac_const = int q^3 F(q) dq = 7/120*pi^4
    real(dl), parameter  :: const2 = 5._dl/7._dl/const_pi**2   !0.072372274_dl
    real(dl), parameter  :: zeta3  = 1.2020569031595942853997_dl
    real(dl), parameter  :: zeta5  = 1.0369277551433699263313_dl
    real(dl), parameter  :: zeta7  = 1.0083492773819228268397_dl

    ! zeta3*3/2/pi^2*4/11*((k_B*COBE_CMBTemp/hbar/c)^3* 8*pi*G/3/(100*km/s/megaparsec)^2/(c^2/eV)
    real(dl), parameter :: neutrino_mass_fac= 94.07_dl !converts omnuh2 into sum m_nu in eV

    integer, parameter  :: nrhopn=2000
    real(dl), parameter :: am_min = 0.01_dl  !0.02_dl
    !smallest a*m_nu to integrate distribution function rather than using series
    real(dl), parameter :: am_max = 600._dl
    !max a*m_nu to integrate

    real(dl), parameter  :: am_minp=am_min*1.1
    real(dl), parameter :: am_maxp=am_max*0.9

    Type TNuPerturbations
        !Sample for massive neutrino momentum
        !Default settings appear to be OK for P_k accuate at 1e-3 level
        integer nqmax !actual number of q modes evolves
        real(dl), allocatable ::  nu_q(:), nu_int_kernel(:)
    contains
    procedure :: init => TNuPerturbations_init
    end type TNuPerturbations

    Type TThermalNuBackground
        !Quantities for the neutrino background momentum distribution assuming thermal
        real(dl) dlnam
        real(dl), dimension(:), allocatable ::  r1,p1,dr1,dp1,ddr1
    contains
    procedure :: init => ThermalNuBackground_init
    procedure :: rho_P => ThermalNuBackground_rho_P
    procedure :: rho => ThermalNuBackground_rho
    procedure :: drho => ThermalNuBackground_drho
    end type TThermalNuBackground

    Type(TThermalNuBackground), save :: ThermalNuBackground

    public fermi_dirac_const,  zeta3, sum_mnu_for_m1, neutrino_mass_fac, TNuPerturbations, ThermalNuBackground
    contains
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine sum_mnu_for_m1(summnu,dsummnu, m1, targ, sgn)
    use constants
    real(dl), intent(in) :: m1, targ, sgn
    real(dl), intent(out) :: summnu, dsummnu
    real(dl) :: m2,m3

    m2 = sqrt(m1**2 + delta_mnu21)
    m3 = sqrt(m1**2 + sgn*delta_mnu31)
    summnu = m1 + m2 + m3 - targ
    dsummnu = m1/m2+m1/m3 + 1

    end subroutine sum_mnu_for_m1

    subroutine TNuPerturbations_init(this,Accuracy)
    !Set up which momenta to integrate the neutrino perturbations, depending on accuracy
    !Using three optimized momenta works very well in most cases
    class(TNuPerturbations) :: this
    real(dl), intent(in) :: Accuracy
    real(dl) :: dq,dlfdlq, q
    integer i

    this%nqmax=3
    if (Accuracy>1) this%nqmax=4
    if (Accuracy>2) this%nqmax=5
    if (Accuracy>3) this%nqmax=nint(Accuracy*10)
    !note this may well be worse than the 5 optimized points

    !We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
    !Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
    !see CAMB notes
    if (allocated(this%nu_q)) deallocate(this%nu_q, this%nu_int_kernel)
    allocate(this%nu_q(this%nqmax))
    allocate(this%nu_int_kernel(this%nqmax))

    if (this%nqmax==3) then
        !Accurate at 2e-4 level
        this%nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
        this%nu_int_kernel(1:3) = (/0.0687359, 3.31435, 2.29911/)
    else if (this%nqmax==4) then
        !This seems to be very accurate (limited by other numerics)
        this%nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
        this%nu_int_kernel(1:4) = (/0.0200251, 1.84539, 3.52736, 0.289427/)
    else if (this%nqmax==5) then
        !exact for n=-4,-2..3
        !This seems to be very accurate (limited by other numerics)
        this%nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)
        this%nu_int_kernel(1:5) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/)
    else
        dq = (12 + this%nqmax/5)/real(this%nqmax)
        do i=1,this%nqmax
            q=(i-0.5d0)*dq
            this%nu_q(i) = q
            dlfdlq=-q/(1._dl+exp(-q))
            this%nu_int_kernel(i)=dq*q**3/(exp(q)+1._dl) * (-0.25_dl*dlfdlq) !now evolve 4F_l/dlfdlq(i)
        end do
    end if
    this%nu_int_kernel=this%nu_int_kernel/fermi_dirac_const

    end subroutine TNuPerturbations_init

    subroutine ThermalNuBackground_init(this)
    class(TThermalNuBackground) :: this
    !  Initialize interpolation tables for massive neutrinos.
    !  Use cubic splines interpolation of log rhonu and pnu vs. log a*m.
    integer i
    real(dl) am, rhonu,pnu
    real(dl) spline_data(nrhopn)

    if (allocated(this%r1)) return

    allocate(this%r1(nrhopn),this%p1(nrhopn),this%dr1(nrhopn),this%dp1(nrhopn),this%ddr1(nrhopn))
    this%dlnam=-(log(am_min/am_max))/(nrhopn-1)

    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), &
    !$OMP& PRIVATE(am,rhonu,pnu)
    do i=1,nrhopn
        am=am_min*exp((i-1)*this%dlnam)
        call nuRhoPres(am,rhonu,pnu)
        this%r1(i)=log(rhonu)
        this%p1(i)=log(pnu)
    end do
    !$OMP END PARALLEL DO

    call splini(spline_data,nrhopn)
    call splder(this%r1,this%dr1,nrhopn,spline_data)
    call splder(this%p1,this%dp1,nrhopn,spline_data)
    call splder(this%dr1,this%ddr1,nrhopn,spline_data)

    end subroutine ThermalNuBackground_init

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine nuRhoPres(am,rhonu,pnu)
    !  Compute the density and pressure of one eigenstate of massive neutrinos,
    !  in units of the mean density of one flavor of massless neutrinos.
    real(dl),  parameter :: qmax=30._dl
    integer, parameter :: nq=100
    real(dl) dum1(nq+1),dum2(nq+1)
    real(dl), intent(in) :: am
    real(dl), intent(out) ::  rhonu,pnu
    integer i
    real(dl) q,aq,v,aqdn,adq

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    !  Integrate up to qmax and then use asymptotic expansion for remainder.
    adq=qmax/nq
    dum1(1)=0._dl
    dum2(1)=0._dl
    do  i=1,nq
        q=i*adq
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        aqdn=adq*q*q*q/(exp(q)+1._dl)
        dum1(i+1)=aqdn/v
        dum2(i+1)=aqdn*v
    end do
    call splint(dum1,rhonu,nq+1)
    call splint(dum2,pnu,nq+1)
    !  Apply asymptotic corrrection for q>qmax and normalize by relativistic
    !  energy density.
    rhonu=(rhonu+dum1(nq+1)/adq)/fermi_dirac_const
    pnu=(pnu+dum2(nq+1)/adq)/fermi_dirac_const/3._dl

    end subroutine nuRhoPres

    !cccccccccccccccccccccccccccccccccccccccccc
    subroutine ThermalNuBackground_rho_P(this,am,rhonu,pnu)
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu, pnu
    real(dl) d
    integer i
    !  Compute massive neutrino density and pressure in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use cubic splines to
    !  interpolate from a table.

    if (am <= am_minp) then
        rhonu=1._dl + const2*am**2
        pnu=(2-rhonu)/3._dl
        return
    else if (am >= am_maxp) then
        rhonu = 3/(2*fermi_dirac_const)*(zeta3*am + (15*zeta5)/2/am)
        pnu = 900._dl/120._dl/fermi_dirac_const*(zeta5-63._dl/4*Zeta7/am**2)/am
        return
    end if

    d=log(am/am_min)/this%dlnam+1._dl
    i=int(d)
    d=d-i

    !  Cubic spline interpolation.
    rhonu=this%r1(i)+d*(this%dr1(i)+d*(3._dl*(this%r1(i+1)-this%r1(i))-2._dl*this%dr1(i) &
        -this%dr1(i+1)+d*(this%dr1(i)+this%dr1(i+1)+2._dl*(this%r1(i)-this%r1(i+1)))))
    pnu=this%p1(i)+d*(this%dp1(i)+d*(3._dl*(this%p1(i+1)-this%p1(i))-2._dl*this%dp1(i) &
        -this%dp1(i+1)+d*(this%dp1(i)+this%dp1(i+1)+2._dl*(this%p1(i)-this%p1(i+1)))))
    rhonu=exp(rhonu)
    pnu=exp(pnu)

    end subroutine ThermalNuBackground_rho_P

    !cccccccccccccccccccccccccccccccccccccccccc
    subroutine ThermalNuBackground_rho(this,am,rhonu)
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu
    real(dl) d
    integer i

    !  Compute massive neutrino density in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use cubic splines to
    !  interpolate from a table.

    if (am <= am_minp) then
        rhonu=1._dl + const2*am**2
        return
    else if (am >= am_maxp) then
        rhonu = 3/(2*fermi_dirac_const)*(zeta3*am + (15*zeta5)/2/am)
        return
    end if

    d=log(am/am_min)/this%dlnam+1._dl
    i=int(d)
    d=d-i

    !  Cubic spline interpolation.
    rhonu=this%r1(i)+d*(this%dr1(i)+d*(3._dl*(this%r1(i+1)-this%r1(i))-2._dl*this%dr1(i) &
        -this%dr1(i+1)+d*(this%dr1(i)+this%dr1(i+1)+2._dl*(this%r1(i)-this%r1(i+1)))))
    rhonu=exp(rhonu)
    end subroutine ThermalNuBackground_rho

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    function ThermalNuBackground_drho(this,am,adotoa,rhonu) result (rhonudot)
    class(TThermalNuBackground) :: this

    !  Compute the time derivative of the mean density in massive neutrinos
    !  and the shear perturbation.
    real(dl) adotoa,rhonu,rhonudot
    real(dl) d
    real(dl), intent(IN) :: am
    integer i

    if (am< am_minp) then
        rhonudot = 2*const2*am**2*adotoa
    else if (am>am_maxp) then
        rhonudot = 3/(2*fermi_dirac_const)*(zeta3*am - (15*zeta5)/2/am)*adotoa
    else
        d=log(am/am_min)/this%dlnam+1._dl
        i=int(d)
        d=d-i
        !  Cubic spline interpolation for rhonudot.
        rhonudot=this%dr1(i)+d*(this%ddr1(i)+d*(3._dl*(this%dr1(i+1)-this%dr1(i)) &
            -2._dl*this%ddr1(i)-this%ddr1(i+1)+d*(this%ddr1(i)+this%ddr1(i+1) &
            +2._dl*(this%dr1(i)-this%dr1(i+1)))))

        rhonudot=rhonu*adotoa*rhonudot/this%dlnam
    end if

    end function ThermalNuBackground_drho

    end module MassiveNu
