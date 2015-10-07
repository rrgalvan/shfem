module nrtypes
  implicit none
  integer, parameter :: i1b = selected_int_kind(2)  ! at least -10^2...10^20
  integer, parameter :: short = i1b
  integer, parameter :: i2b = selected_int_kind(4)  ! at least -10^4...10^4
  integer, parameter :: med = i2b
  integer, parameter :: i4b = selected_int_kind(9)  ! at least -10^9...10^9
  integer, parameter :: long = i4b

  integer, parameter :: sp = kind(1.0) ! single precision
  integer, parameter :: dp = kind(1.0d0) ! double precission

  real(sp), parameter :: pi=3.141592653589793238462643383279502884197_sp
  real(sp), parameter :: pio2=1.57079632679489661923132169163975144209858_sp
  real(sp), parameter :: twopi=6.283185307179586476925286766559005768394_sp
  real(sp), parameter :: sqrt2=1.41421356237309504880168872420969807856967_sp
  real(sp), parameter :: euler=0.5772156649015328606065120900824024310422_sp
  real(dp), parameter :: pi_dp=3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: pio2_dp=1.57079632679489661923132169163975144209858_dp
  real(dp), parameter :: twopi_dp=6.283185307179586476925286766559005768394_dp
end module nrtypes
