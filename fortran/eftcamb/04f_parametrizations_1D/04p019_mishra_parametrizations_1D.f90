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

!> @file 04p019_mishra_parametrizations_1D.f90
!! This file contains the definition of the Mishra parametrization, around a=0,
!! up to third order, inheriting from parametrized_function_1D.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the Mishra parametrization, around a=0,
!! up to third order, inheriting from parametrized_function_1D.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_mishra_parametrizations_1D

    use precision
    use EFT_def
    use MpiUtils
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public mishra_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the Mishra parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: mishra_parametrization_1D

        real(dl) :: Omega_V
        real(dl) :: lambda0

    contains

        ! utility functions:
        procedure :: set_param_number      => MishraParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the Mishra parametrized function.
        procedure :: init_parameters       => MishraParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => MishraParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => MishraParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => MishraParameterNames                    !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => MishraParameterNamesLatex               !< subroutine that returns the i-th parameter name of the function in latex format.

        ! evaluation procedures:
        procedure :: value                 => MishraParametrized1DValue               !< function that returns the value of the Mishra.
        procedure :: first_derivative      => MishraParametrized1DFirstDerivative     !< function that returns the first derivative of the Mishra.
        procedure :: second_derivative     => MishraParametrized1DSecondDerivative    !< function that returns the second derivative of the Mishra.
        procedure :: third_derivative      => MishraParametrized1DThirdDerivative     !< function that returns the third derivative of the Mishra.
        procedure :: integral              => MishraParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type mishra_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the Mishra parametrization.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the Mishra parametrized function.
    subroutine MishraParametrized1DSetParamNumber( self )

        implicit none

        class(mishra_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 2

    end subroutine MishraParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine MishraParametrized1DInitParams( self, array )

        implicit none

        class(mishra_parametrization_1D)                        :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%Omega_V = array(1)
        self%lambda0 = array(2)

    end subroutine MishraParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine MishraParametrized1DParameterValues( self, i, value )

        implicit none

        class(mishra_parametrization_1D)       :: self        !< the base class
        integer     , intent(in)               :: i           !< The index of the parameter
        real(dl)    , intent(out)              :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%Omega_V
            case(2)
                value = self%lambda0
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine MishraParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine MishraParametrized1DFeedback( self, print_params )

        implicit none

        class(mishra_parametrization_1D) :: self         !< the base class
        logical, optional                :: print_params !< optional flag that decised whether to print numerical values
                                                         !! of the parameters.

        integer                                 :: i
        real(dl)                                :: param_value
        character(len=EFT_names_max_length)     :: param_name
        logical                                 :: print_params_temp

        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        write(*,*)     'Mishra parametrization: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine MishraParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine MishraParameterNames( self, i, name )

        implicit none

        class(mishra_parametrization_1D)  :: self   !< the base class
        integer     , intent(in)           :: i      !< the index of the parameter
        character(*), intent(out)          :: name   !< the output name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'In parametrized_function_1D:', self%name
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the parameter name:
        if ( allocated(self%param_names) ) then
            name = self%param_names(i)%string
        else
            if ( i==1 ) then
                name = TRIM(self%name)//'_omega_v'
            else if ( i==2 ) then
                name = TRIM(self%name)//'_lambda0'
            end if
        end if

    end subroutine MishraParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine MishraParameterNamesLatex( self, i, latexname )

        implicit none

        class(mishra_parametrization_1D)   :: self       !< the base class
        integer     , intent(in)            :: i          !< the index of the parameter
        character(*), intent(out)           :: latexname  !< the output latex name of the i-th parameter

        ! check the input index:
        if ( i>self%parameter_number ) then
            write(*,*) 'In parametrized_function_1D:', self%name
            write(*,*) 'Illegal index for parameter_names.'
            write(*,*) 'Maximum value is:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if
        ! return the parameter name:
        if ( allocated(self%param_names_latex) ) then
            latexname = self%param_names_latex(i)%string
        else
            if ( i==1 ) then
                latexname = TRIM(self%name_latex)//'_{\omega_v}'
            else if ( i==2 ) then
                latexname = TRIM(self%name_latex)//'_{\lambda_0}'
            end if
        end if

    end subroutine MishraParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function MishraParametrized1DValue( self, x, eft_cache )

        implicit none

        class(mishra_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: MishraParametrized1DValue                           !< the output value
        real(dl) :: A0
        real(dl) :: y0
        real(dl) :: ZZ
        A0 = 1._dl/self%Omega_V-1._dl
        ZZ = x**3._dl+A0
        y0 = (x**3._dl/ZZ)**0.5_dl

        MishraParametrized1DValue = -(y0 + self%lambda0/3._dl/ZZ*self%Omega_V**(-0.5_dl)*(-y0*A0*Log(self%Omega_V**(-0.5_dl)-1._dl) + -x**3._dl + y0/self%Omega_V**0.5_dl + y0*A0*Log(ZZ**0.5_dl-x**1.5_dl)))**2._dl

    end function MishraParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function MishraParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(mishra_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: MishraParametrized1DFirstDerivative                 !< the output value
        real(dl) :: A0
        real(dl) :: y0
        real(dl) :: AA
        real(dl) :: BB
        real(dl) :: yy
        real(dl) :: yy1
        real(dl) :: ZZ
        A0 = 1._dl/self%Omega_V-1._dl
        ZZ = x**3._dl+A0
        y0 = (x**3._dl/ZZ)**0.5_dl
        AA = A0-2._dl*x**3._dl
        BB = (1._dl+A0*x**(-3._dl))**(-0.5_dl)
        yy = y0 + self%lambda0/3._dl/ZZ*self%Omega_V**(-0.5_dl)*(-y0*A0*Log(self%Omega_V**(-0.5_dl)-1._dl) + -x**3._dl + y0/self%Omega_V**0.5_dl + y0*A0*Log(ZZ**0.5_dl-x**1.5_dl))
        yy1 = 1.5_dl*A0*x**(-3._dl)*y0**3._dl + 0.5_dl*self%Omega_V**(-0.5_dl)*ZZ**(-2._dl)*self%lambda0*(-A0*y0*AA*Log(self%Omega_V**(-0.5_dl)-1._dl) + self%Omega_V**(-0.5_dl)*(y0*(AA-3*A0*self%Omega_V**0.5_dl*BB*ZZ) + A0*self%Omega_V**0.5_dl*BB*AA*Log(ZZ**0.5_dl-x**1.5_dl)))

        MishraParametrized1DFirstDerivative = -2/x*yy*yy1

    end function MishraParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function MishraParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(mishra_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: MishraParametrized1DSecondDerivative                !< the output value
        real(dl) :: A0
        real(dl) :: y0
        real(dl) :: AA
        real(dl) :: BB
        real(dl) :: yy
        real(dl) :: yy1
        real(dl) :: ZZ
        real(dl) :: yy2
        A0 = 1._dl/self%Omega_V-1._dl
        ZZ = x**3._dl+A0
        y0 = (x**3._dl/ZZ)**0.5_dl
        AA = A0-2._dl*x**3._dl
        BB = (1._dl+A0*x**(-3._dl))**(-0.5_dl)
        yy = y0 + self%lambda0/3._dl/ZZ*self%Omega_V**(-0.5_dl)*(-y0*A0*Log(self%Omega_V**(-0.5_dl)-1._dl) + -x**3._dl + y0/self%Omega_V**0.5_dl + y0*A0*Log(ZZ**0.5_dl-x**1.5_dl))
        yy1 = 1.5_dl*A0*x**(-3._dl)*y0**3._dl + 0.5_dl*self%Omega_V**(-0.5_dl)*ZZ**(-2._dl)*self%lambda0*(-A0*y0*AA*Log(self%Omega_V**(-0.5_dl)-1._dl) + self%Omega_V**(-0.5_dl)*(y0*(AA-3*A0*self%Omega_V**0.5_dl*BB*ZZ) + A0*self%Omega_V**0.5_dl*BB*AA*Log(ZZ**0.5_dl-x**1.5_dl)))
        yy2 = (3._dl * y0)/(4._dl * ZZ**3) * (4._dl * x**6 * (1._dl + A0) * self%lambda0 + A0**2 * (3._dl * ZZ + self%lambda0 * (1._dl + A0 - 7._dl * (1._dl + A0)**0.5_dl * y0 * ZZ)) - 2._dl * x**3 * A0 * (3._dl * ZZ + self%lambda0 * (5._dl + 5._dl * A0 - 4._dl * (1._dl + A0)**0.5_dl * y0 * ZZ)) - A0 * (1._dl + A0)**0.5_dl * (4._dl * x**6 - 10._dl * x**3 * A0 + A0**2) * self%lambda0 * (log(-1._dl + (1._dl + A0)**0.5_dl) - log(-x**(3._dl/2._dl) + ZZ**0.5_dl)))

        MishraParametrized1DSecondDerivative = 2._dl*x**(-2._dl)*(yy*yy1 - yy1**2._dl - yy*yy2)

    end function MishraParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function MishraParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(mishra_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: MishraParametrized1DThirdDerivative                 !< the output value
        real(dl) :: A0
        real(dl) :: y0
        real(dl) :: AA
        real(dl) :: BB
        real(dl) :: yy
        real(dl) :: yy1
        real(dl) :: ZZ
        real(dl) :: yy2
        real(dl) :: yy3a
        A0 = 1._dl/self%Omega_V-1._dl
        ZZ = x**3._dl+A0
        y0 = (x**3._dl/ZZ)**0.5_dl
        AA = A0-2._dl*x**3._dl
        BB = (1._dl+A0*x**(-3._dl))**(-0.5_dl)
        yy = y0 + self%lambda0/3._dl/ZZ*self%Omega_V**(-0.5_dl)*(-y0*A0*Log(self%Omega_V**(-0.5_dl)-1._dl) + -x**3._dl + y0/self%Omega_V**0.5_dl + y0*A0*Log(ZZ**0.5_dl-x**1.5_dl))
        yy1 = 1.5_dl*A0*x**(-3._dl)*y0**3._dl + 0.5_dl*self%Omega_V**(-0.5_dl)*ZZ**(-2._dl)*self%lambda0*(-A0*y0*AA*Log(self%Omega_V**(-0.5_dl)-1._dl) + self%Omega_V**(-0.5_dl)*(y0*(AA-3*A0*self%Omega_V**0.5_dl*BB*ZZ) + A0*self%Omega_V**0.5_dl*BB*AA*Log(ZZ**0.5_dl-x**1.5_dl)))
        yy2 = (3._dl * y0)/(4._dl * ZZ**3) * (4._dl * x**6 * (1._dl + A0) * self%lambda0 + A0**2 * (3._dl * ZZ + self%lambda0 * (1._dl + A0 - 7._dl * (1._dl + A0)**0.5_dl * y0 * ZZ)) - 2._dl * x**3 * A0 * (3._dl * ZZ + self%lambda0 * (5._dl + 5._dl * A0 - 4._dl * (1._dl + A0)**0.5_dl * y0 * ZZ)) - A0 * (1._dl + A0)**0.5_dl * (4._dl * x**6 - 10._dl * x**3 * A0 + A0**2) * self%lambda0 * (log(-1._dl + (1._dl + A0)**0.5_dl) - log(-x**(3._dl/2._dl) + ZZ**0.5_dl)))
        yy3a = 9._dl/(8._dl * ZZ**4) * (-4._dl * x**9 * self%lambda0 * (5._dl * A0 * (1._dl + A0)**0.5_dl + 2._dl * y0 + 2._dl * A0 * y0) + A0**3 * y0 * (self%lambda0 + A0 * self%lambda0 + 3._dl * ZZ) + 2._dl * x**6 * A0 * (35._dl * A0 * (1._dl + A0)**0.5_dl * self%lambda0 + 30._dl * (1._dl + A0) * self%lambda0 * y0 + 6._dl * y0 * ZZ) - 3._dl * x**3 * A0**2 * (5._dl * A0 * (1._dl + A0)**0.5_dl * self%lambda0 + 12._dl * (1._dl + A0) * self%lambda0 * y0 + 10._dl * y0 * ZZ) - A0 * (1._dl + A0)**0.5_dl * (-8._dl * x**9 + 60._dl * x**6 * A0 - 36._dl * x**3 * A0**2 + A0**3) * self%lambda0 * y0 * (log(-1._dl + (1._dl + A0)**0.5_dl) - log(-x**(3._dl/2._dl) + ZZ**0.5_dl)))

        MishraParametrized1DThirdDerivative = -2._dl*x**(-3._dl)*(3._dl*yy1*(yy2-yy1) + yy*(yy3a+2._dl*yy1-3._dl*yy2))

    end function MishraParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function MishraParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(mishra_parametrization_1D)                   :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: MishraParametrized1DIntegral                        !< the output value
        real(dl) :: A0
        real(dl) :: y0
        real(dl) :: ZZ
        A0 = 1._dl/self%Omega_V-1._dl
        ZZ = x**3._dl+A0
        y0 = (x**3._dl/ZZ)**0.5_dl

        MishraParametrized1DIntegral = (x**3._dl + A0)/x/(A0+1)*Exp(-2._dl/3._dl*self%Omega_V**(-0.5_dl)*self%lambda0/ZZ*(self%Omega_V**(-0.5_dl)-(x**3._dl*(A0+x**3._dl))**0.5_dl+(A0+2._dl*x**3._dl)*(Log(self%Omega_V**(-0.5_dl)-1)-ZZ*Log(ZZ**0.5_dl-x**1.5_dl))))

    end function MishraParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_mishra_parametrizations_1D

!----------------------------------------------------------------------------------------
