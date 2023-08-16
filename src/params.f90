module params

  use iso_fortran_env

  implicit none
  
  !number of threads when executing in parallel or using OpenBLAS
  integer, parameter :: max_num_threads=20
  
  !real precision parameter
  integer, parameter :: rp=real64
  !integer precision parameter
  integer, parameter :: ip=int32

  !unit number for the program log and unit offset for other units
  !use unit_log=6 to output to terminal, otherwise make unit_log greater than 10
  !offset_units should be greater than unit_log
  integer, parameter :: unit_log=6, offset_units=20

  real(rp), parameter :: pi=4.0d0*datan(1.0d0)
  real(rp), parameter :: two_pi=8.0d0*datan(1.0d0)
  real(rp), parameter :: half_pi=2.0d0*datan(1.0d0)
  real(rp), parameter :: third_pi=4.0d0*datan(1.0d0)/3.0d0
  real(rp), parameter :: quarter_pi=datan(1.0d0)

contains

end module params
