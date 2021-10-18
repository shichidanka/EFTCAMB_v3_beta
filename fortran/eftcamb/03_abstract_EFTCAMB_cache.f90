!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2020 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 03_abstract_EFTCAMB_cache.f90
!! This file contains the definition of the EFTCAMB caches.
!! These are used to store parameters that can be used by EFTCAMB, in TEFTCAMB_parameter_cache
!! or are used to store partial results when solving the time evolution of perturbations,
!! in TEFTCAMB_timestep_cache.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the EFTCAMB caches.
!! These are used to store parameters that can be used by EFTCAMB, in TEFTCAMB_parameter_cache
!! or are used to store partial results when solving the time evolution of perturbations,
!! in TEFTCAMB_timestep_cache.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_cache

    use precision
    use MassiveNu
    use classes
    use MpiUtils
    use FileUtils

    implicit none

    private

    public TEFTCAMB_parameter_cache, TEFTCAMB_timestep_cache

    ! some settings:
    character(*), parameter :: cache_output_format = 'es19.10e3'

    !define debug file
    type(TTextFile) :: file_FRW
    type(TTextFile) :: file_BDens
    type(TTextFile) :: file_BPres
    type(TTextFile) :: file_BOmegas
    type(TTextFile) :: file_BackgroundEFT
    type(TTextFile) :: file_SecondOrdEFT
    type(TTextFile) :: file_BackgroundQ
    type(TTextFile) :: file_EinsteinCoeff
    type(TTextFile) :: file_PiCoeff
    type(TTextFile) :: file_PiSolution
    type(TTextFile) :: file_EinsteinSol
    type(TTextFile) :: file_TensorCoeff
    type(TTextFile) :: file_Sources
    type(TTextFile) :: file_MetricMG

    !----------------------------------------------------------------------------------------
    !> This is the type that defines the EFTCAMB parameter cache. The idea is to copy in here
    !! the cosmological parameters that we need from CAMB and then use this for the interfaces.
    !! This also contains a wrapper for Massive Neutrinos quantities that can be used in
    !! the background code.
    type, extends(TCambComponent) :: TEFTCAMB_parameter_cache

        ! 1) relative densities:
        real(dl) :: omegac         !< the value of \f$ \Omega_{\rm CDM}^0 \f$.
        real(dl) :: omegab         !< the value of \f$ \Omega_{\rm b}^0 \f$.
        real(dl) :: omegav         !< the value of \f$ \Omega_{\Lambda}^0 \f$.
        real(dl) :: omegak         !< the value of \f$ \Omega_{\rm K}^0 \f$.
        real(dl) :: omegan         !< the value of \f$ \Omega_{m\nu}^0 \f$ of massive neutrinos.
        real(dl) :: omegag         !< the value of \f$ \Omega_{ \gamma }^0 \f$.
        real(dl) :: omegar         !< the value of \f$ \Omega_{ \nu }^0 \f$ of massless neutrinos.
        ! 2) Hubble constant:
        real(dl) :: h0             !< reduced Hubble constant \f$ H_0/100 \f$.
        real(dl) :: h0_Mpc         !< the Hubble constant in MegaParsec \f$ 10^3 \cdot H_0/c \f$.
        ! 3) densities:
        real(dl) :: grhog          !< the value of \f$ 8 \pi G_{N} \rho_{\gamma}(t_0) \f$.
        real(dl) :: grhornomass    !< the value of \f$ 8 \pi G_{N} \rho_{\nu}(t_0) \f$.
        real(dl) :: grhoc          !< the value of \f$ 8 \pi G_{N} \rho_{\rm CDM}(t_0) \f$.
        real(dl) :: grhob          !< the value of \f$ 8 \pi G_{N} \rho_{\rm b}(t_0) \f$.
        real(dl) :: grhov          !< the value of \f$ 8 \pi G_{N} \rho_{\Lambda}(t_0) \f$.
        real(dl) :: grhok          !< the value of \f$ 8 \pi G_{N} \rho_{\rm K}(t_0) \f$.
        ! 4) massive neutrinos:
        integer  :: Num_Nu_Massive                        !< number of massive neutrinos
        integer  :: Nu_mass_eigenstates                   !< number of mass eigenstates
        real(dl), allocatable, dimension(:) :: grhormass  !< densities of neutrinos in each mass eigenstate \f$ 8 \pi G_{N} \rho_{ m\nu }(t_0) \f$
        real(dl), allocatable, dimension(:) :: nu_masses  !< neutrino masses

    contains

        procedure :: initialize => EFTCAMBParameterCacheInit  !< subroutine that initializes to zero all the elements of the parameter cache.
        procedure :: is_nan     => EFTCAMBParameterCacheIsNan !< Subroutine that check if an element of the TEFTCAMB_parameter_cache is Nan.
        procedure :: print      => EFTCAMBParameterCachePrint !< subroutine that prints the EFTCAMB parameters cache to screen.
        ! python interface:
        procedure, nopass :: SelfPointer => TEFTCAMB_parameter_cache_SelfPointer

    end type TEFTCAMB_parameter_cache

    !----------------------------------------------------------------------------------------
    !> This is the type that defines the EFTCAMB time step cache.
    type, extends(TCambComponent) :: TEFTCAMB_timestep_cache

        ! 1) time and k:
        real(dl) :: a             !< the value of the scale factor at which the cache is being used.
        real(dl) :: tau           !< the value of conformal time at the given scale factor.
        real(dl) :: k             !< the scale that is being solved for. In \f$ \mbox{Mpc}^{-1} \f$.
        ! 2) total matter densities:
        real(dl) :: grhoa2        !< the input value of \f$ \sum_m\rho_m / a^2 m_0^2 \f$.
        real(dl) :: grhom_t       !< the value of \f$ \sum_m\rho_m a^2 /m_0^2 \f$.
        real(dl) :: gpresm_t      !< the value of \f$ \sum_m P_m a^2 /m_0^2 \f$.
        real(dl) :: gpresdotm_t   !< the value of \f$ \sum_m\dot{P}_m a^2 /m_0^2 \f$.
        ! 3) densities and pressure of the various species:
        real(dl) :: grhob_t       !< the value of \f$ \rho_b a^2 / m_0^2 \f$
        real(dl) :: grhoc_t       !< the value of \f$ \rho_{cdm} a^2 / m_0^2 \f$
        real(dl) :: grhor_t       !< the value of \f$ \rho_{\nu} a^2 / m_0^2 \f$
        real(dl) :: grhog_t       !< the value of \f$ \rho_{\gamma} a^2 / m_0^2 \f$
        real(dl) :: grhov_t       !< the value of \f$ \rho_{\Lambda} a^2 / m_0^2 \f$. Used if needed, especially by designer models.
        real(dl) :: gpiv_t        !< the value of \f$ \P_{\Lambda} a^2 / m_0^2 \f$. Used if needed, especially by designer models.
        real(dl) :: grhonu_tot    !< the value of \f$ \sum_\nu \rho_{m\nu} a^2 / m_0^2 \f$
        real(dl) :: gpinu_tot     !< the value of \f$ \sum_\nu P_{m\nu} a^2 / m_0^2 \f$
        real(dl) :: grhonudot_tot !< the value of \f$ \sum_\nu \dot{\rho}_{m\nu} a^2 / m_0^2 \f$
        real(dl) :: gpinudot_tot  !< the value of \f$ \sum_\nu \dot{P}_{m\nu} a^2 / m_0^2 \f$
        ! 4) expansion history:
        real(dl) :: adotoa        !< the value of \f$ \mathcal{H} \f$ at the given scale factor.
        real(dl) :: Hdot          !< the value of \f$ d\mathcal{H} /d \tau \f$ at the given scale factor.
        real(dl) :: Hdotdot       !< the value of \f$ d^2 \mathcal{H} / d \tau^2 \f$ at the given scale factor.
        ! 5) EFT functions:
        real(dl) :: EFTOmegaV     !< the value of Omega \f$ \Omega(a) \f$.
        real(dl) :: EFTOmegaP     !< the value of the derivative wrt scale factor of Omega \f$ d \Omega(a) / da \f$.
        real(dl) :: EFTOmegaPP    !< the value of the second derivative wrt scale factor of Omega \f$ d^2 \Omega(a) / da^2 \f$.
        real(dl) :: EFTOmegaPPP   !< the value of the third derivative wrt scale factor of Omega \f$ d^3 \Omega(a) / da^3 \f$.
        real(dl) :: EFTc          !< the value of \f$ c a^2/m_0^2 \f$.
        real(dl) :: EFTcdot       !< the value of \f$ \dot{c} a^2/m_0^2 \f$.
        real(dl) :: EFTcdotdot    !< the value of \f$ \dot{c} a^2/m_0^2 \f$.
        real(dl) :: EFTLambda     !< the value of \f$ \Lambda a^2/m_0^2 \f$.
        real(dl) :: EFTLambdadot  !< the value of \f$ \dot{\Lambda}a^2/m_0^2 \f$. Derivative of \f$ \Lambda\f$ wrt conformal time.
        real(dl) :: EFTGamma1V    !< the value of Gamma 1 \f$ \gamma_1(a) \f$.
        real(dl) :: EFTGamma1P    !< the value of the derivative wrt scale factor of Gamma 1 \f$  d \gamma_1(a) / da \f$.
        real(dl) :: EFTGamma1PP   !< the value of the derivative wrt scale factor of Gamma 1 \f$  d \gamma_1(a) / da \f$.
        real(dl) :: EFTGamma2V    !< the value of Gamma 2 \f$ \gamma_2(a) \f$.
        real(dl) :: EFTGamma2P    !< the value of the derivative wrt scale factor of Gamma 2 \f$  d \gamma_2(a) / da \f$.
        real(dl) :: EFTGamma2PP   !< the value of the derivative wrt scale factor of Gamma 2 \f$  d \gamma_2(a) / da \f$.
        real(dl) :: EFTGamma3V    !< the value of Gamma 3 \f$ \gamma_3(a) \f$.
        real(dl) :: EFTGamma3P    !< the value of the derivative wrt scale factor of Gamma 3 \f$  d \gamma_3(a) / da \f$.
        real(dl) :: EFTGamma3PP   !< the value of the derivative wrt scale factor of Gamma 3 \f$  d \gamma_3(a) / da \f$.
        real(dl) :: EFTGamma4V    !< the value of Gamma 4 \f$ \gamma_4(a) \f$.
        real(dl) :: EFTGamma4P    !< the value of the derivative wrt scale factor of Gamma 4 \f$  d \gamma_4(a) / da \f$.
        real(dl) :: EFTGamma4PP   !< the value of the second derivative wrt scale factor of Gamma 4 \f$  d^2 \gamma_4(a) / da^2 \f$.
        real(dl) :: EFTGamma5V    !< the value of Gamma 5 \f$ \gamma_5(a) \f$.
        real(dl) :: EFTGamma5P    !< the value of the derivative wrt scale factor of Gamma 5 \f$  d \gamma_5(a) / da \f$.
        real(dl) :: EFTGamma6V    !< the value of Gamma 6 \f$ \gamma_6(a) \f$.
        real(dl) :: EFTGamma6P    !< the value of the derivative wrt scale factor of Gamma 6 \f$  d \gamma_6(a) / da \f$.
        ! 6) other background quantities:
        real(dl) :: grhoq         !< the value of the effective density of the Q field. Refer to the Numerical Notes for the definition.
        real(dl) :: gpresq        !< the value of the effective pressure of the Q field. Refer to the Numerical Notes for the definition.
        real(dl) :: grhodotq      !< the value of the time derivative of the effective density of the Q field. Refer to the Numerical Notes for the definition.
        real(dl) :: gpresdotq     !< the value of the time derivative of the effective pressure of the Q field. Refer to the Numerical Notes for the definition.
        ! 7) the Einstein equations coefficients:
        real(dl) :: EFTeomF       !< the value of the Einstein equations coefficient F. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomN       !< the value of the Einstein equations coefficient N. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomNdot    !< the value of the Einstein equations coefficient dN/dtau. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomX       !< the value of the Einstein equations coefficient X. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomXdot    !< the value of the Einstein equations coefficient dX/dtau. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomY       !< the value of the Einstein equations coefficient Y. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomG       !< the value of the Einstein equations coefficient G. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomU       !< the value of the Einstein equations coefficient U. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomL       !< the value of the Einstein equations coefficient L. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomM       !< the value of the Einstein equations coefficient M. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomV       !< the value of the Einstein equations coefficient V. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomVdot    !< the value of the Einstein equations coefficient dV/dtau. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTeomQ       !< the value of the Einstein equations coefficient Q. Refer to the Numerical Notes for the definition.
        ! 8) pi field factors:
        real(dl) :: EFTpiA1       !< the value of the pi field equation coefficient A1. Scale independent part of A. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiA2       !< the value of the pi field equation coefficient A2. Part proportional to \f$ k^2 \f$. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiB1       !< the value of the pi field equation coefficient B1. Scale independent part of B. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiB2       !< the value of the pi field equation coefficient B2. Part proportional to \f$ k^2 \f$. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiC        !< the value of the pi field equation coefficient C. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiD1       !< the value of the pi field equation coefficient D1. Scale independent part of D. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiD2       !< the value of the pi field equation coefficient D2. Part proportional to \f$ k^2 \f$. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTpiE        !< the value of the pi field equation coefficient E. Refer to the Numerical Notes for the definition.
        ! 9) pi field quantities:
        real(dl) :: pi            !< the value of the pi field at a given time and scale.
        real(dl) :: pidot         !< the value of the (conformal) time derivative of the pi field at a given time and scale.
        real(dl) :: pidotdot      !< the value of the (conformal) second time derivative of the pi field at a given time and scale.
        ! 10) scalar perturbations quantities:
        real(dl) :: etak          !< Syncronous gauge \f$ eta*k \f$ perturbation.
        real(dl) :: etakdot       !< Syncronous gauge \f$ \dot{eta}*k \f$ perturbation. This is the Einstein equation that is actually integrated.
        real(dl) :: z             !< Syncronous gauge Z perturbation.
        real(dl) :: dz            !< Syncronous gauge dot Z perturbation.
        real(dl) :: sigma         !< Syncronous gauge sigma perturbation.
        real(dl) :: sigmadot      !< Syncronous gauge dot sigma perturbation.
        real(dl) :: sigmadotdot   !< Syncronous gauge second time derivative of sigma perturbation.
        real(dl) :: clxc          !< Syncronous gauge cdm density perturbation.
        real(dl) :: clxb          !< Syncronous gauge baryon density perturbation.
        real(dl) :: clxg          !< Syncronous gauge radiation density perturbation.
        real(dl) :: clxr          !< Syncronous gauge massless neutrinos density perturbation.
        real(dl) :: vb            !< Syncronous gauge baryon velocity.
        real(dl) :: dgpnu         !< Syncronous gauge massive neutrinos pressure perturbation.
        real(dl) :: dgrho         !< Syncronous gauge total matter density perturbation.
        real(dl) :: dgq           !< Syncronous gauge total matter velocity perturbation.
        ! 11) tensor perturbations quantities:
        real(dl) :: EFTAT         !< the value of the tensor equation coefficient A. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTBT         !< the value of the tensor equation coefficient B. Refer to the Numerical Notes for the definition.
        real(dl) :: EFTDT         !< the value of the tensor equation coefficient D. Refer to the Numerical Notes for the definition.
        ! 12) Kinetic and Gradient quantities for the stability check:
        real(dl) :: EFT_kinetic   !< the value of the kinetic term. Refer to the Numerical Notes for the definition.
        real(dl) :: EFT_gradient  !< the value of the gradient term. Refer to the Numerical Notes for the definition.
        real(dl) :: EFT_mu1       !< the value of the first mass eigenvalue. Refer to the Numerical Notes for the definition.
        real(dl) :: EFT_mu2       !< the value of the second mass eigenvalue. Refer to the Numerical Notes for the definition.
        ! 13) other quantities usefull for debug purposes:
        real(dl) :: Psi           !< Perturbation in the 00 component of the metric in Conformal Newtonian gauge.
        real(dl) :: Phi           !< Perturbation in the space component of the metric in Conformal Newtonian gauge.
        real(dl) :: PsiDot        !< Time derivatives of the Newtonian gauge metric perturbations.
        real(dl) :: PhiDot        !< Time derivatives of the Newtonian gauge metric perturbations.
        real(dl) :: EFTISW        !< Source for ISW effect.
        real(dl) :: EFTLensing    !< Source for lensing effect.
        real(dl) :: T_monopole    !< The effective CMB temperature monopole.
        real(dl) :: T_dipole      !< The effective CMB temperature dipole.
        real(dl) :: mu            !< Effective perturbation gravitational constant. \f$ k^2\Psi/\Delta_{m} \f$
        real(dl) :: gamma         !< Ratio between the two gravitational potentials. \f$ \gamma = \Phi/\Psi \f$
        real(dl) :: sigma_eft     !< Ratio between Weyl potential and matter density. \f$ \Sigma = \frac{\Phi+\Psi}{\Delta_{m}} \f$
        real(dl) :: dgrho_Newt    !< Newtonian gauge total matter density perturbation.
        real(dl) :: delta_DE_eff  !< Syncronous gauge Dark energy fluid perturbation.

    contains

        procedure :: initialize        => EFTCAMBTimestepCacheInit      !< subroutine that initializes to zero all the elements of the cache.
        procedure :: is_nan            => EFTCAMBTimestepCacheIsNan     !< Subroutine that check if an element of the TEFTCAMB_timestep_cache is Nan.
        procedure :: open_cache_files  => EFTCAMBTimestepCacheOpenFile  !< subroutine that opens the files to dump the cache to file.
        procedure :: dump_cache_files  => EFTCAMBTimestepCacheDumpFile  !< subroutine that dumps the cache to the files.
        procedure :: close_cache_files => EFTCAMBTimestepCacheCloseFile !< subroutine that closes the files where the cache has benn dumped.
        ! python interface:
        procedure, nopass :: SelfPointer => TEFTCAMB_timestep_cache_SelfPointer

    end type TEFTCAMB_timestep_cache

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes to zero all the elements of the cache.
    subroutine EFTCAMBTimestepCacheInit( self )

        implicit none

        class(TEFTCAMB_timestep_cache)  :: self !< the base class.

        ! initialize all class members to zero:
        ! 1) time and k:
        self%a             = 0._dl
        self%tau           = 0._dl
        self%k             = 0._dl
        ! 2) expansion history:
        self%adotoa        = 0._dl
        self%Hdot          = 0._dl
        self%Hdotdot       = 0._dl
        ! 3) total matter densities:
        self%grhoa2        = 0._dl
        self%grhom_t       = 0._dl
        self%gpresm_t      = 0._dl
        self%gpresdotm_t   = 0._dl
        ! 4) densities and pressure of the various species:
        self%grhob_t       = 0._dl
        self%grhoc_t       = 0._dl
        self%grhor_t       = 0._dl
        self%grhog_t       = 0._dl
        self%grhov_t       = 0._dl
        self%gpiv_t        = 0._dl
        self%grhonu_tot    = 0._dl
        self%gpinu_tot     = 0._dl
        self%grhonudot_tot = 0._dl
        self%gpinudot_tot  = 0._dl
        ! 5) EFT functions:
        self%EFTOmegaV     = 0._dl
        self%EFTOmegaP     = 0._dl
        self%EFTOmegaPP    = 0._dl
        self%EFTOmegaPPP   = 0._dl
        self%EFTc          = 0._dl
        self%EFTcdot       = 0._dl
        self%EFTLambda     = 0._dl
        self%EFTLambdadot  = 0._dl
        self%EFTGamma1V    = 0._dl
        self%EFTGamma1P    = 0._dl
        self%EFTGamma2V    = 0._dl
        self%EFTGamma2P    = 0._dl
        self%EFTGamma3V    = 0._dl
        self%EFTGamma3P    = 0._dl
        self%EFTGamma4V    = 0._dl
        self%EFTGamma4P    = 0._dl
        self%EFTGamma4PP   = 0._dl
        self%EFTGamma5V    = 0._dl
        self%EFTGamma5P    = 0._dl
        self%EFTGamma6V    = 0._dl
        self%EFTGamma6P    = 0._dl
        ! 6) other background quantities:
        self%grhoq         = 0._dl
        self%gpresq        = 0._dl
        self%grhodotq      = 0._dl
        self%gpresdotq     = 0._dl
        ! 7) the Einstein equations coefficients:
        self%EFTeomF       = 0._dl
        self%EFTeomN       = 0._dl
        self%EFTeomNdot    = 0._dl
        self%EFTeomX       = 0._dl
        self%EFTeomXdot    = 0._dl
        self%EFTeomY       = 0._dl
        self%EFTeomG       = 0._dl
        self%EFTeomU       = 0._dl
        self%EFTeomL       = 0._dl
        self%EFTeomM       = 0._dl
        self%EFTeomV       = 0._dl
        self%EFTeomVdot    = 0._dl
        self%EFTeomQ       = 0._dl
        ! 8) pi field factors:
        self%EFTpiA1       = 0._dl
        self%EFTpiA2       = 0._dl
        self%EFTpiB1       = 0._dl
        self%EFTpiB2       = 0._dl
        self%EFTpiC        = 0._dl
        self%EFTpiD1       = 0._dl
        self%EFTpiD2       = 0._dl
        self%EFTpiE        = 0._dl
        ! 9) pi field quantities:
        self%pi            = 0._dl
        self%pidot         = 0._dl
        self%pidotdot      = 0._dl
        ! 10) perturbations quantities:
        self%etak          = 0._dl
        self%etakdot       = 0._dl
        self%z             = 0._dl
        self%dz            = 0._dl
        self%sigma         = 0._dl
        self%sigmadot      = 0._dl
        self%sigmadotdot   = 0._dl
        self%clxc          = 0._dl
        self%clxb          = 0._dl
        self%clxg          = 0._dl
        self%clxr          = 0._dl
        self%vb            = 0._dl
        self%dgpnu         = 0._dl
        self%dgrho         = 0._dl
        self%dgq           = 0._dl
        ! 11) tensor perturbations quantities:
        self%EFTAT         = 0._dl
        self%EFTBT         = 0._dl
        self%EFTDT         = 0._dl
        ! 12) Kinetic and Gradient quantities for the stability check:
        self%EFT_kinetic   = 0._dl
        self%EFT_gradient  = 0._dl
        self%EFT_mu1       = 0._dl
        self%EFT_mu2       = 0._dl
        ! 13) other quantities usefull for debug purposes:
        self%Psi           = 0._dl
        self%Phi           = 0._dl
        self%PsiDot        = 0._dl
        self%PhiDot        = 0._dl
        self%EFTISW        = 0._dl
        self%EFTLensing    = 0._dl
        self%T_monopole    = 0._dl
        self%T_dipole      = 0._dl
        self%mu            = 0._dl
        self%gamma         = 0._dl
        self%dgrho_Newt    = 0._dl
        self%delta_DE_eff  = 0._dl

    end subroutine EFTCAMBTimestepCacheInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that check if an element of the TEFTCAMB_timestep_cache is Nan.
    subroutine EFTCAMBTimestepCacheIsNan( self, HaveNan )

        implicit none

        class(TEFTCAMB_timestep_cache), intent(in)  :: self    !< The base class.
        logical, intent(out)                       :: HaveNan !< Logical variable which describes the presence of a Nan variable.
                                                              !< If an element of the TEFTCAMB_timestep_cache is Nan, you get HaveNan=.True.

        HaveNan = .False.
        HaveNan = HaveNan.or.IsNaN(self%a)
        HaveNan = HaveNan.or.IsNaN(self%tau)
        HaveNan = HaveNan.or.IsNaN(self%k)
        HaveNan = HaveNan.or.IsNaN(self%grhoa2)
        HaveNan = HaveNan.or.IsNaN(self%grhom_t)
        HaveNan = HaveNan.or.IsNaN(self%gpresm_t)
        HaveNan = HaveNan.or.IsNaN(self%gpresdotm_t)
        HaveNan = HaveNan.or.IsNaN(self%grhob_t)
        HaveNan = HaveNan.or.IsNaN(self%grhoc_t)
        HaveNan = HaveNan.or.IsNaN(self%grhor_t)
        HaveNan = HaveNan.or.IsNaN(self%grhog_t)
        HaveNan = HaveNan.or.IsNaN(self%grhov_t)
        HaveNan = HaveNan.or.IsNaN(self%gpiv_t)
        HaveNan = HaveNan.or.IsNaN(self%grhonu_tot)
        HaveNan = HaveNan.or.IsNaN(self%gpinu_tot)
        HaveNan = HaveNan.or.IsNaN(self%grhonudot_tot)
        HaveNan = HaveNan.or.IsNaN(self%gpinudot_tot)
        HaveNan = HaveNan.or.IsNaN(self%adotoa)
        HaveNan = HaveNan.or.IsNaN(self%Hdot)
        HaveNan = HaveNan.or.IsNaN(self%Hdotdot)
        HaveNan = HaveNan.or.IsNaN(self%EFTOmegaV)
        HaveNan = HaveNan.or.IsNaN(self%EFTOmegaP)
        HaveNan = HaveNan.or.IsNaN(self%EFTOmegaPP)
        HaveNan = HaveNan.or.IsNaN(self%EFTOmegaPPP)
        HaveNan = HaveNan.or.IsNaN(self%EFTc)
        HaveNan = HaveNan.or.IsNaN(self%EFTcdot)
        HaveNan = HaveNan.or.IsNaN(self%EFTcdotdot)
        HaveNan = HaveNan.or.IsNaN(self%EFTLambda)
        HaveNan = HaveNan.or.IsNaN(self%EFTLambdadot)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma1V)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma1P)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma1PP)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma2V)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma2P)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma2PP)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma3V)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma3P)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma3PP)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma4V)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma4P)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma4PP)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma5V)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma5P)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma6V)
        HaveNan = HaveNan.or.IsNaN(self%EFTGamma6P)
        HaveNan = HaveNan.or.IsNaN(self%grhoq)
        HaveNan = HaveNan.or.IsNaN(self%gpresq)
        HaveNan = HaveNan.or.IsNaN(self%grhodotq)
        HaveNan = HaveNan.or.IsNaN(self%gpresdotq)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomF)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomN)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomNdot)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomX)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomXdot)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomY)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomG)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomU)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomL)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomM)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomV)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomVdot)
        HaveNan = HaveNan.or.IsNaN(self%EFTeomQ)
        HaveNan = HaveNan.or.IsNaN(self%EFTpiA1)
        HaveNan = HaveNan.or.IsNaN(self%EFTpiA2)
        HaveNan = HaveNan.or.IsNaN(self%EFTpiB1)
        HaveNan = HaveNan.or.IsNaN(self%EFTpiB2)
        HaveNan = HaveNan.or.IsNaN(self%EFTpiC)
        HaveNan = HaveNan.or.IsNaN(self%EFTpiD1)
        HaveNan = HaveNan.or.IsNaN(self%EFTpiD2)
        HaveNan = HaveNan.or.IsNaN(self%EFTpiE)
        HaveNan = HaveNan.or.IsNaN(self%pi)
        HaveNan = HaveNan.or.IsNaN(self%pidot)
        HaveNan = HaveNan.or.IsNaN(self%pidotdot)
        HaveNan = HaveNan.or.IsNaN(self%etak)
        HaveNan = HaveNan.or.IsNaN(self%etakdot)
        HaveNan = HaveNan.or.IsNaN(self%z)
        HaveNan = HaveNan.or.IsNaN(self%dz)
        HaveNan = HaveNan.or.IsNaN(self%sigma)
        HaveNan = HaveNan.or.IsNaN(self%sigmadot)
        HaveNan = HaveNan.or.IsNaN(self%sigmadotdot)
        HaveNan = HaveNan.or.IsNaN(self%clxc)
        HaveNan = HaveNan.or.IsNaN(self%clxb)
        HaveNan = HaveNan.or.IsNaN(self%clxg)
        HaveNan = HaveNan.or.IsNaN(self%clxr)
        HaveNan = HaveNan.or.IsNaN(self%vb)
        HaveNan = HaveNan.or.IsNaN(self%dgpnu)
        HaveNan = HaveNan.or.IsNaN(self%dgrho)
        HaveNan = HaveNan.or.IsNaN(self%dgq)
        HaveNan = HaveNan.or.IsNaN(self%EFTAT)
        HaveNan = HaveNan.or.IsNaN(self%EFTBT)
        HaveNan = HaveNan.or.IsNaN(self%EFTDT)
        HaveNan = HaveNan.or.IsNaN(self%EFT_kinetic)
        HaveNan = HaveNan.or.IsNaN(self%EFT_gradient)
        HaveNan = HaveNan.or.IsNaN(self%Psi)
        HaveNan = HaveNan.or.IsNaN(self%Phi)
        HaveNan = HaveNan.or.IsNaN(self%PsiDot)
        HaveNan = HaveNan.or.IsNaN(self%PhiDot)
        HaveNan = HaveNan.or.IsNaN(self%EFTISW)
        HaveNan = HaveNan.or.IsNaN(self%EFTLensing)
        HaveNan = HaveNan.or.IsNaN(self%T_monopole)
        HaveNan = HaveNan.or.IsNaN(self%T_dipole)
        HaveNan = HaveNan.or.IsNaN(self%mu)
        HaveNan = HaveNan.or.IsNaN(self%gamma)
        HaveNan = HaveNan.or.IsNaN(self%dgrho_Newt)
        HaveNan = HaveNan.or.IsNaN(self%delta_DE_eff)

    end subroutine EFTCAMBTimestepCacheIsNan

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that opens the files to dump the cache to file.
    subroutine EFTCAMBTimestepCacheOpenFile( self, outroot )

        implicit none

        class(TEFTCAMB_timestep_cache)  :: self    !< the base class.
        character(len=*), intent(in)   :: outroot !< output root of the file.

        logical :: is_open

        ! print some feedback:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') 'EFTCAMB cache opening print files in:'
        write(*,'(a)')  outroot
        write(*,'(a)') "***************************************************************"

        !set output units
        file_FRW%unit           = 111
        file_BDens%unit         = 222
        file_BPres%unit         = 333
        file_BOmegas%unit       = 444
        file_BackgroundEFT%unit = 555
        file_SecondOrdEFT%unit  = 666
        file_BackgroundQ%unit   = 777
        file_EinsteinCoeff%unit = 888
        file_PiCoeff%unit       = 999
        file_PiSolution%unit    = 1111
        file_EinsteinSol%unit   = 2222
        file_TensorCoeff%unit   = 3333
        file_Sources%unit       = 4444
        file_MetricMG%unit      = 5555

        ! test whether the files are already open:
        call test_open( file_FRW%unit           )
        call test_open( file_BDens%unit         )
        call test_open( file_BPres%unit         )
        call test_open( file_BOmegas%unit       )
        call test_open( file_BackgroundEFT%unit )
        call test_open( file_SecondOrdEFT%unit  )
        call test_open( file_BackgroundQ%unit   )
        call test_open( file_EinsteinCoeff%unit )
        call test_open( file_PiCoeff%unit       )
        call test_open( file_PiSolution%unit    )
        call test_open( file_EinsteinSol%unit   )
        call test_open( file_TensorCoeff%unit   )
        call test_open( file_Sources%unit       )
        call test_open( file_MetricMG%unit      )

        ! open the files:
        call file_FRW%CreateFile(            TRIM(outroot)//'cache_FRW.dat'            )
        call file_BDens%CreateFile(          TRIM(outroot)//'cache_BDens.dat'          )
        call file_BPres%CreateFile(          TRIM(outroot)//'cache_BPres.dat'          )
        call file_BOmegas%CreateFile(        TRIM(outroot)//'cache_BOmegas.dat'        )
        call file_BackgroundEFT%CreateFile(  TRIM(outroot)//'cache_BackgroundEFT.dat'  )
        call file_SecondOrdEFT%CreateFile(   TRIM(outroot)//'cache_SecondOrdEFT.dat'   )
        call file_BackgroundQ%CreateFile(    TRIM(outroot)//'cache_BackgroundQ.dat'    )
        call file_EinsteinCoeff%CreateFile(  TRIM(outroot)//'cache_EinsteinCoeff.dat'  )
        call file_PiCoeff%CreateFile(        TRIM(outroot)//'cache_PiCoeff.dat'        )
        call file_PiSolution%CreateFile(     TRIM(outroot)//'cache_PiSolution.dat'     )
        call file_EinsteinSol%CreateFile(    TRIM(outroot)//'cache_EinsteinSol.dat'    )
        call file_TensorCoeff%CreateFile(    TRIM(outroot)//'cache_TensorCoeff.dat'    )
        call file_Sources%CreateFile(        TRIM(outroot)//'cache_Sources.dat'        )
        call file_MetricMG%CreateFile(       TRIM(outroot)//'cache_MetricMG.dat'       )

        ! write the headers:
        write (file_FRW%unit          ,'(12a)') '# ', 'a ', 'tau ', 'k ', 'adotoa ', 'Hdot ', 'Hdotdot '
        write (file_BDens%unit        ,'(12a)') '# ', 'a ', 'tau ', 'k ', 'grhom_t ', 'grhob_t ', 'grhoc_t ', 'grhor_t ', 'grhog_t ', 'grhov_t ', 'grhonu_tot ', 'grhonudot_tot '
        write (file_BPres%unit        ,'(12a)') '# ', 'a ', 'tau ', 'k ', 'gpresm_t ', 'gpresdotm_t ', 'gpiv_t ', 'gpinu_tot ', 'gpinudot_tot '
        write (file_BOmegas%unit      ,'(12a)') '# ', 'a ', 'tau ', 'k ', 'omegam ', 'omegab ', 'omegac ', 'omegar ', 'omegag ', 'omegav ', 'omeganu_tot '
        write (file_BackgroundEFT%unit,'(20a)') '# ', 'a ', 'tau ', 'k ', 'EFTOmegaV ', 'EFTOmegaP ', 'EFTOmegaPP ', 'EFTOmegaPPP ', 'EFTc ', 'EFTcdot ', 'EFTLambda ', 'EFTLambdadot '
        write (file_SecondOrdEFT%unit ,'(20a)') '# ', 'a ', 'tau ', 'k ', 'EFTGamma1V ', 'EFTGamma1P ', 'EFTGamma2V ', 'EFTGamma2P ', 'EFTGamma3V ', 'EFTGamma3P ', 'EFTGamma4V ', 'EFTGamma4P ', 'EFTGamma4PP ', 'EFTGamma5V ', 'EFTGamma5P ', 'EFTGamma6V ', 'EFTGamma6P '
        write (file_BackgroundQ%unit  ,'(12a)') '# ', 'a ', 'tau ', 'k ', 'grhoq ', 'gpresq ', 'grhodotq ', 'gpresdotq '
        write (file_EinsteinCoeff%unit,'(18a)') '# ', 'a ', 'tau ', 'k ', 'EFTeomF ', 'EFTeomN ', 'EFTeomNdot ', 'EFTeomX ', 'EFTeomXdot ', 'EFTeomY ', 'EFTeomG ', 'EFTeomU ', 'EFTeomL ', 'EFTeomM ', 'EFTeomV ', 'EFTeomVdot ', 'EFTeomQ'
        write (file_PiCoeff%unit      ,'(14a)') '# ', 'a ', 'tau ', 'k ', 'EFTpiA1 ', 'EFTpiA2 ', 'EFTpiB1 ', 'EFTpiB2 ', 'EFTpiC ', 'EFTpiD1 ', 'EFTpiD2 ', 'EFTpiE '
        write (file_PiSolution%unit   ,'(12a)') '# ', 'a ', 'tau ', 'k ', 'pi ', 'pidot ', 'pidotdot '
        write (file_EinsteinSol%unit  ,'(20a)') '# ', 'a ', 'tau ', 'k ', 'z ', 'sigma ', 'clxc ', 'clxb ', 'clxg ', 'clxr ', 'dgpnu ', 'dgrho ', 'dgq '
        write (file_TensorCoeff%unit  ,'(12a)') '# ', 'a ', 'tau ', 'k ', 'EFTAT ', 'EFTBT ', 'EFTDT '
        write (file_Sources%unit      ,'(12a)') '# ', 'a ', 'tau ', 'k ', 'EFTISW ', 'EFTLensing ', 'Theta0 ', 'Theta1 '
        write (file_MetricMG%unit     ,'(12a)') '# ', 'a ', 'tau ', 'k ', 'Psi ', 'Phi ', 'PsiDot ', 'PhiDot ', 'mu ', 'gamma ', 'delta_DE_eff '

    contains

        ! Temporary subroutine that tests wether a cache file is open.
        subroutine test_open( number )

            implicit none

            integer, intent(in) :: number
            logical             :: is_open

            inquire( unit=number, opened=is_open )
            if ( is_open ) then
                write(*,*) 'EFTCAMB ERROR: Oputput unit', number, 'is already open.'
                write(*,*) 'EFTCAMB cannot use it and cannot proceed.'
                call MpiStop('EFTCAMB error')
            end if

        end subroutine test_open

    end subroutine EFTCAMBTimestepCacheOpenFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that closes the files where the cache has benn dumped.
    subroutine EFTCAMBTimestepCacheCloseFile( self )

        implicit none

        class(TEFTCAMB_timestep_cache)  :: self !< the base class.

        ! test wether the files can be closed:
        call test_close( file_FRW%unit          )
        call test_close( file_BDens%unit        )
        call test_close( file_BPres%unit        )
        call test_close( file_BOmegas%unit      )
        call test_close( file_BackgroundEFT%unit)
        call test_close( file_SecondOrdEFT%unit )
        call test_close( file_BackgroundQ%unit  )
        call test_close( file_EinsteinCoeff%unit)
        call test_close( file_PiCoeff%unit      )
        call test_close( file_PiSolution%unit   )
        call test_close( file_EinsteinSol%unit  )
        call test_close( file_TensorCoeff%unit  )
        call test_close( file_Sources%unit      )
        call test_close( file_MetricMG%unit     )

        ! close the files:
        call file_FRW%close()
        call file_BDens%close()
        call file_BPres%close()
        call file_BOmegas%close()
        call file_BackgroundEFT%close()
        call file_SecondOrdEFT%close()
        call file_BackgroundQ%close()
        call file_EinsteinCoeff%close()
        call file_PiCoeff%close()
        call file_PiSolution%close()
        call file_EinsteinSol%close()
        call file_TensorCoeff%close()
        call file_Sources%close()
        call file_MetricMG%close()

        ! print some feedback:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') 'EFTCAMB cache printing done.'
        write(*,'(a)') "***************************************************************"

    contains

        ! Temporary subroutine that tests wether a cahce file is open.
        subroutine test_close( number )

            implicit none

            integer, intent(in) :: number
            logical             :: is_open

            inquire( unit=number, opened=is_open )
            if ( .not. is_open ) then
                write(*,*) 'EFTCAMB ERROR: Oputput unit', number, 'is not open.'
                write(*,*) 'EFTCAMB is trying to close it and cannot proceed.'
                call MpiStop('EFTCAMB error')
            end if

        end subroutine test_close

    end subroutine EFTCAMBTimestepCacheCloseFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that dumps the cache to the files.
    subroutine EFTCAMBTimestepCacheDumpFile( self )

        implicit none

        class(TEFTCAMB_timestep_cache)  :: self !< the base class.

        ! write the background expansion history:
        write (file_FRW%unit           ,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%adotoa, self%Hdot, self%Hdotdot
        ! write the background densities:
        write (file_BDens%unit         ,'(14'//cache_output_format//')')  self%a, self%tau, self%k, self%grhom_t, self%grhob_t, self%grhoc_t, self%grhor_t, self%grhog_t, self%grhov_t, self%grhonu_tot, self%grhonudot_tot
        ! write the background pressure:
        write (file_BPres%unit         ,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%gpresm_t, self%gpresdotm_t, self%gpiv_t, self%gpinu_tot, self%gpinudot_tot
        ! write the background omegas:
        write (file_BOmegas%unit       ,'(14'//cache_output_format//')')  self%a, self%tau, self%k, self%grhom_t/(3._dl*self%adotoa**2), self%grhob_t/(3._dl*self%adotoa**2), self%grhoc_t/(3._dl*self%adotoa**2), self%grhor_t/(3._dl*self%adotoa**2), self%grhog_t/(3._dl*self%adotoa**2), self%grhov_t/(3._dl*self%adotoa**2), self%grhonu_tot/(3._dl*self%adotoa**2)
        ! write background EFT functions:
        write (file_BackgroundEFT%unit ,'(14'//cache_output_format//')')  self%a, self%tau, self%k, self%EFTOmegaV, self%EFTOmegaP, self%EFTOmegaPP, self%EFTOmegaPPP, self%EFTc, self%EFTcdot, self%EFTLambda, self%EFTLambdadot
        ! write second order EFT functions:
        write (file_SecondOrdEFT%unit  ,'(18'//cache_output_format//')')  self%a, self%tau, self%k, self%EFTGamma1V, self%EFTGamma1P, self%EFTGamma2V, self%EFTGamma2P, self%EFTGamma3V, self%EFTGamma3P, self%EFTGamma4V, self%EFTGamma4P, self%EFTGamma4PP, self%EFTGamma5V, self%EFTGamma5P, self%EFTGamma6V, self%EFTGamma6P
        ! write background EFT auxiliary quantities:
        write (file_BackgroundQ%unit   ,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%grhoq, self%gpresq, self%grhodotq, self%gpresdotq
        ! write Einstein equations coefficients:
        write (file_EinsteinCoeff%unit ,'(18'//cache_output_format//')')  self%a, self%tau, self%k, self%EFTeomF, self%EFTeomN, self%EFTeomNdot, self%EFTeomX, self%EFTeomXdot, self%EFTeomY, self%EFTeomG, self%EFTeomU, self%EFTeomL, self%EFTeomM, self%EFTeomV, self%EFTeomVdot, self%EFTeomQ
        ! write pi field coefficients:
        write (file_PiCoeff%unit       ,'(14'//cache_output_format//')')  self%a, self%tau, self%k, self%EFTpiA1, self%EFTpiA2, self%EFTpiB1, self%EFTpiB2, self%EFTpiC, self%EFTpiD1, self%EFTpiD2, self%EFTpiE
        ! write pi field solution:
        write (file_PiSolution%unit    ,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%pi, self%pidot, self%pidotdot
        ! write some perturbations:
        write (file_EinsteinSol%unit   ,'(20'//cache_output_format//')')  self%a, self%tau, self%k, self%z, self%sigma, self%clxc, self%clxb, self%clxg, self%clxr, self%dgpnu, self%dgrho, self%dgq
        ! write tensor coefficients:
        write (file_TensorCoeff%unit   ,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%EFTAT, self%EFTBT, self%EFTDT
        ! write sources:
        write (file_Sources%unit       ,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%EFTISW, self%EFTLensing, self%T_monopole, self%T_dipole
        ! write metric potentials:
        write (file_MetricMG%unit      ,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%Psi, self%Phi, self%PhiDot, self%PsiDot, self%mu, self%gamma, self%delta_DE_eff

    end subroutine EFTCAMBTimestepCacheDumpFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine for the python interface:
    subroutine TEFTCAMB_timestep_cache_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TEFTCAMB_timestep_cache), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P
        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TEFTCAMB_timestep_cache_SelfPointer

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes to zero all the elements of the parameter cache.
    subroutine EFTCAMBParameterCacheInit( self )

        implicit none

        class(TEFTCAMB_parameter_cache)  :: self !< the base class.

        ! initialize all class members to zero:
        self%omegac      = 0._dl
        self%omegab      = 0._dl
        self%omegav      = 0._dl
        self%omegak      = 0._dl
        self%omegan      = 0._dl
        self%omegag      = 0._dl
        self%omegar      = 0._dl
        self%h0          = 0._dl
        self%h0_Mpc      = 0._dl
        self%grhog       = 0._dl
        self%grhornomass = 0._dl
        self%grhoc       = 0._dl
        self%grhob       = 0._dl
        self%grhov       = 0._dl
        self%grhok       = 0._dl
        self%Num_Nu_Massive       = 0
        self%Nu_mass_eigenstates  = 0

    end subroutine EFTCAMBParameterCacheInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints the EFTCAMB parameters cache to screen.
    subroutine EFTCAMBParameterCachePrint( self )

        implicit none

        class(TEFTCAMB_parameter_cache)  :: self !< the base class.

        integer  :: i

        ! print to screen the parameter cache:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') 'EFTCAMB parameters cache content:'
        write(*,'(a)') "***************************************************************"
        write(*,'(a14,E13.6)') ' Omega_CDM  : ', self%omegac
        write(*,'(a14,E13.6)') ' Omega_b    : ', self%omegab
        write(*,'(a14,E13.6)') ' Omega_v    : ', self%omegav
        write(*,'(a14,E13.6)') ' Omega_k    : ', self%omegak
        write(*,'(a14,E13.6)') ' Omega_n    : ', self%omegan
        write(*,'(a14,E13.6)') ' Omega_g    : ', self%omegag
        write(*,'(a14,E13.6)') ' Omega_r    : ', self%omegar
        write(*,'(a14,F12.6)') ' h          : ', self%h0
        write(*,'(a14,E13.6)') ' h_Mpc      : ', self%h0_Mpc
        write(*,'(a14,E13.6)') ' grhog      : ', self%grhog
        write(*,'(a14,E13.6)') ' grnonomass : ', self%grhornomass
        write(*,'(a14,E13.6)') ' grhoc      : ', self%grhoc
        write(*,'(a14,E13.6)') ' grhob      : ', self%grhob
        write(*,'(a14,E13.6)') ' grhov      : ', self%grhov
        write(*,'(a14,E13.6)') ' grhok      : ', self%grhok
        write(*,'(a22,I10)') ' Num_Nu_Massive      : ', self%Num_Nu_Massive
        write(*,'(a22,I10)') ' Nu_mass_eigenstates : ', self%Nu_mass_eigenstates
        do i=1, self%Nu_mass_eigenstates
            write(*,'(a11,I3,a9,E13.6)') ' grhormass(',i,')      : ', self%grhormass(i)
            write(*,'(a11,I3,a9,E13.6)') ' nu_masses(',i,')      : ', self%nu_masses(i)
        end do
        write(*,'(a)') "***************************************************************"

    end subroutine EFTCAMBParameterCachePrint

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that check if an element of the TEFTCAMB_parameter_cache is Nan.
    subroutine EFTCAMBParameterCacheIsNan( self, HaveNan )

        implicit none

        class(TEFTCAMB_parameter_cache), intent(in)  :: self    !< The base class.
        logical, intent(out)                        :: HaveNan !< Logical variable which describes the presence of a Nan variable.
                                                               !< If an element of the TEFTCAMB_parameter_cache is Nan, you get HaveNan=.True.

        integer                                     :: i

        HaveNan = .False.
        HaveNan = HaveNan.or.IsNaN(self%omegac)
        HaveNan = HaveNan.or.IsNaN(self%omegab)
        HaveNan = HaveNan.or.IsNaN(self%omegav)
        HaveNan = HaveNan.or.IsNaN(self%omegak)
        HaveNan = HaveNan.or.IsNaN(self%omegan)
        HaveNan = HaveNan.or.IsNaN(self%omegag)
        HaveNan = HaveNan.or.IsNaN(self%omegar)
        HaveNan = HaveNan.or.IsNaN(self%h0)
        HaveNan = HaveNan.or.IsNaN(self%h0_Mpc)
        HaveNan = HaveNan.or.IsNaN(self%grhog)
        HaveNan = HaveNan.or.IsNaN(self%grhornomass)
        HaveNan = HaveNan.or.IsNaN(self%grhoc)
        HaveNan = HaveNan.or.IsNaN(self%grhob)
        HaveNan = HaveNan.or.IsNaN(self%grhov)
        HaveNan = HaveNan.or.IsNaN(self%grhok)
        HaveNan = HaveNan.or.IsNaN(self%Num_Nu_Massive*1.0_dl)
        HaveNan = HaveNan.or.IsNaN(self%Nu_mass_eigenstates*1.0)

        do i=1, self%Nu_mass_eigenstates
            HaveNan = HaveNan.or.IsNaN(self%grhormass(i)).or.IsNaN(self%nu_masses(i))
        end do

    end subroutine EFTCAMBParameterCacheIsNan

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine for the python interface:
    subroutine TEFTCAMB_parameter_cache_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TEFTCAMB_parameter_cache), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P
        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TEFTCAMB_parameter_cache_SelfPointer

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_cache

!----------------------------------------------------------------------------------------
