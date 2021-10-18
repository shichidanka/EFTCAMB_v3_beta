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

!> @file 04p018_steplog_parametrizations_1D.f90
!! This file contains the definition of the steplog parametrization,
!! inheriting from parametrized_function_1D.

!----------------------------------------------------------------------------------------
!> This module contains the definition of the steplog parametrization,
!! inheriting from parametrized_function_1D.

!> @author Marco Raveri, Meng-Xiang Lin

module EFTCAMB_steplog_parametrizations_1D

    use precision
    use MpiUtils
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_parametrizations_1D

    implicit none

    private

    public steplog_parametrization_1D

    ! ---------------------------------------------------------------------------------------------
    !> Type containing the steplog function parametrization. Inherits from parametrized_function_1D.
    type, extends ( parametrized_function_1D ) :: steplog_parametrization_1D

        real(dl) :: v1      ! The value of the function when a is -infinity (in cosmology we use a=0)
        real(dl) :: v2      ! The value of the function when a is infinity (in cosmology we use a=1)
        real(dl) :: aT      ! The time of transition
        real(dl) :: Delta   ! The width of transition

    contains

        ! utility functions:
        procedure :: set_param_number      => steplogParametrized1DSetParamNumber      !< subroutine that sets the number of parameters of the steplog parametrized function.
        procedure :: init_parameters       => steplogParametrized1DInitParams          !< subroutine that initializes the function parameters based on the values found in an input array.
        procedure :: parameter_value       => steplogParametrized1DParameterValues     !< subroutine that returns the value of the function i-th parameter.
        procedure :: feedback              => steplogParametrized1DFeedback            !< subroutine that prints to screen the informations about the function.
        procedure :: parameter_names       => steplogParameterNames                    !< subroutine that returns the i-th parameter name of the function.
        procedure :: parameter_names_latex => steplogParameterNamesLatex               !< subroutine that returns the i-th parameter name of the function in latex format.

        ! evaluation procedures:
        procedure :: value                 => steplogParametrized1DValue               !< function that returns the value of the steplog function.
        procedure :: first_derivative      => steplogParametrized1DFirstDerivative     !< function that returns the first derivative of the steplog function.
        procedure :: second_derivative     => steplogParametrized1DSecondDerivative    !< function that returns the second derivative of the steplog function.
        procedure :: third_derivative      => steplogParametrized1DThirdDerivative     !< function that returns the third derivative of the steplog function.
        procedure :: integral              => steplogParametrized1DIntegral            !< function that returns the strange integral that we need for w_DE.

    end type steplog_parametrization_1D

contains

    ! ---------------------------------------------------------------------------------------------
    ! Implementation of the steplog function.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the number of parameters of the steplog parametrized function.
    subroutine steplogParametrized1DSetParamNumber( self )

        implicit none

        class(steplog_parametrization_1D) :: self       !< the base class

        ! initialize the number of parameters:
        self%parameter_number = 4

    end subroutine steplogParametrized1DSetParamNumber

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the function parameters based on the values found in an input array.
    subroutine steplogParametrized1DInitParams( self, array )

        implicit none

        class(steplog_parametrization_1D)                       :: self   !< the base class.
        real(dl), dimension(self%parameter_number), intent(in)  :: array  !< input array with the values of the parameters.

        self%v1 = array(1)
        self%v2 = array(2)
        self%aT = array(3)
        self%Delta = array(4)

    end subroutine steplogParametrized1DInitParams

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the value of the function i-th parameter.
    subroutine steplogParametrized1DParameterValues( self, i, value )

        implicit none

        class(steplog_parametrization_1D)    :: self        !< the base class
        integer     , intent(in)             :: i           !< The index of the parameter
        real(dl)    , intent(out)            :: value       !< the output value of the i-th parameter

        select case (i)
            case(1)
                value = self%v1
            case(2)
                value = self%v2
            case(3)
                value = self%aT
            case(4)
                value = self%Delta
            case default
                write(*,*) 'Illegal index for parameter_names.'
                write(*,*) 'Maximum value is:', self%parameter_number
                call MpiStop('EFTCAMB error')
        end select

    end subroutine steplogParametrized1DParameterValues

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the informations about the function.
    subroutine steplogParametrized1DFeedback( self, print_params )

        implicit none

        class(steplog_parametrization_1D)    :: self         !< the base class
        logical, optional                    :: print_params !< optional flag that decised whether to print numerical values
                                                             !! of the parameters.

        integer                               :: i
        real(dl)                              :: param_value
        character(len=EFT_names_max_length)   :: param_name
        logical                               :: print_params_temp

        if ( present(print_params) ) then
            print_params_temp = print_params
        else
            print_params_temp = .True.
        end if

        write(*,*)     'steplog parametrization: ', self%name
        if ( print_params_temp ) then
            do i=1, self%parameter_number
                call self%parameter_names( i, param_name  )
                call self%parameter_value( i, param_value )
                write(*,'(a23,a,F12.6)') param_name, '=', param_value
            end do
        end if

    end subroutine steplogParametrized1DFeedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the i-th parameter name.
    subroutine steplogParameterNames( self, i, name )

        implicit none

        class(steplog_parametrization_1D)  :: self   !< the base class
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
                name = TRIM(self%name)//'_v1'
            else if ( i==2 ) then
                name = TRIM(self%name)//'_v2'
            else if ( i==3 ) then
                name = TRIM(self%name)//'_at'
            else if ( i==4 ) then
                name = TRIM(self%name)//'_delta'
            end if
        end if

    end subroutine steplogParameterNames

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that returns the latex version of the i-th parameter name.
    subroutine steplogParameterNamesLatex( self, i, latexname )

        implicit none

        class(steplog_parametrization_1D)   :: self       !< the base class
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
                latexname = TRIM(self%name_latex)//'_1'
            else if ( i==2 ) then
                latexname = TRIM(self%name_latex)//'_2'
            else if ( i==3 ) then
                latexname = TRIM(self%name_latex)//'_{a_T}'
            else if ( i==4 ) then
                latexname = TRIM(self%name_latex)//'_{\delta}'
            end if
        end if

    end subroutine steplogParameterNamesLatex

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the function in the scale factor.
    function steplogParametrized1DValue( self, x, eft_cache )

        implicit none

        class(steplog_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: steplogParametrized1DValue                              !< the output value

        if (x==0) then
            steplogParametrized1DValue = self%v1
        else
            steplogParametrized1DValue = 0.5_dl*(self%v2-self%v1)*((log(x)-log(self%aT))/self%Delta) /sqrt(1+((log(x)-log(self%aT))/self%Delta)**2) + 0.5_dl*(self%v1+self%v2)
        endif

    end function steplogParametrized1DValue

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the value of the first derivative, wrt scale factor, of the function.
    function steplogParametrized1DFirstDerivative( self, x, eft_cache )

        implicit none

        class(steplog_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: steplogParametrized1DFirstDerivative                    !< the output value

        if (x==0) then
            steplogParametrized1DFirstDerivative = 0._dl
        else
            steplogParametrized1DFirstDerivative = 0.5_dl*(self%v2-self%v1)/self%Delta/x /(1+((log(x)-log(self%aT))/self%Delta)**2)**1.5_dl
        endif

    end function steplogParametrized1DFirstDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the second derivative of the function.
    function steplogParametrized1DSecondDerivative( self, x, eft_cache )

        implicit none

        class(steplog_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: steplogParametrized1DSecondDerivative                   !< the output value

        if (x==0) then
            steplogParametrized1DSecondDerivative = 0._dl
        else
            steplogParametrized1DSecondDerivative = 0.5_dl*(self%v1-self%v2)*self%Delta *(self%Delta**2+log(x)**2+log(x)*(3-2*log(self%aT))+(-3+log(self%aT))*log(self%aT)) /x**2 /sqrt(1+((log(x)-log(self%aT))/self%Delta)**2) /((log(x)-log(self%aT))**2+self%Delta**2)**2
        endif

    end function steplogParametrized1DSecondDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the third derivative of the function.
    function steplogParametrized1DThirdDerivative( self, x, eft_cache )

        implicit none

        class(steplog_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: steplogParametrized1DThirdDerivative                    !< the output value
        real(dl) :: term1

        if (x==0) then
            steplogParametrized1DThirdDerivative = 0._dl
        else
            term1 = 9*self%Delta**2 + 2*log(x)**3 + log(x)**2*(9-6*log(self%aT)) - 4*log(self%aT)*(3+self%Delta**2) + 9*log(self%aT)**2 - 2*log(self%aT)**3 + 2*log(x)*(6+2*self%Delta**2+3*log(self%aT)*(-3+log(self%aT)))
            steplogParametrized1DThirdDerivative = 0.5_dl*(self%v2-self%v1)*self%Delta *(self%Delta**2*(-3+2*self%Delta**2)+(log(x)-log(self%aT))*term1) /x**3 /sqrt(1+((log(x)-log(self%aT))/self%Delta)**2) /((log(x)-log(self%aT))**2+self%Delta**2)**3
        endif

    end function steplogParametrized1DThirdDerivative

    ! ---------------------------------------------------------------------------------------------
    !> Function that returns the integral of the function, as defined in the notes.
    function steplogParametrized1DIntegral( self, x, eft_cache )

        implicit none

        class(steplog_parametrization_1D)                      :: self      !< the base class
        real(dl), intent(in)                               :: x         !< the input scale factor
        type(TEFTCAMB_timestep_cache), intent(in), optional :: eft_cache !< the optional input EFTCAMB cache
        real(dl) :: steplogParametrized1DIntegral                           !< the output value
        real(dl) :: term1, term2

        if (x==0) then
            steplogParametrized1DIntegral = 0._dl
        else
            term1 = x**( -1._dl -1.5_dl*(self%v1+self%v2) - 1.5_dl*(self%v2-self%v1)*(log(x)-2*log(self%aT)+log(self%aT)**2/log(x))/sqrt(self%Delta**2+(log(x)-log(self%aT))**2) )
            term2 = exp( 1.5_dl*(self%v1-self%v2) *( self%Delta**2/sqrt(self%Delta**2+(log(x)-log(self%aT))**2) - sqrt(self%Delta**2+log(self%aT)**2) ) )
            steplogParametrized1DIntegral = term1*term2
        endif

    end function steplogParametrized1DIntegral

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_steplog_parametrizations_1D

!----------------------------------------------------------------------------------------
