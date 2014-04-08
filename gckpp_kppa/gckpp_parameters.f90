!---------------------- BEGIN gckpp_parameters.f90 BEGIN ---------------------
! @file gckpp_parameters.f90                                                  
! @author yifanwu                                                             
! @date 2014-04-08 15:56:45.013319                                            
! @brief Program parameters                                                   
!                                                                             
! Integration tolerances, program constants, and species indices.             
!                                                                             
! This file was generated by Kppa: http://www.paratools.com/Kppa              
!-----------------------------------------------------------------------------


MODULE gckpp_parameters


  IMPLICIT NONE


!-----------------------------------------------------------------------------
! Integration tolerances                                                      
!-----------------------------------------------------------------------------

  ! Absolute tolerance 
  REAL(8), PARAMETER :: ATOLS = 1.0
  ! Relative tolerance 
  REAL(8), PARAMETER :: RTOLS = 0.001


!-----------------------------------------------------------------------------
! Concentration constants                                                     
!-----------------------------------------------------------------------------

  ! Conversion factor 
  REAL(8), PARAMETER :: CFACTOR = 1.0
  ! Default initialization value for variable species 
  REAL(8), PARAMETER :: VAR_DEFAULT = 0.0
  ! Default initialization value for fixed species 
  REAL(8), PARAMETER :: FIX_DEFAULT = 0.0


!-----------------------------------------------------------------------------
! Program constants                                                           
!-----------------------------------------------------------------------------

  ! Species count 
  INTEGER, PARAMETER :: NSPEC = 107
  ! Variable species count 
  INTEGER, PARAMETER :: NVAR = 92
  ! Fixed species count 
  INTEGER, PARAMETER :: NFIX = 15
  ! Active variable species count 
  INTEGER, PARAMETER :: NVARACT = 85
  ! Reaction (equation) count 
  INTEGER, PARAMETER :: NREACT = 331


!-----------------------------------------------------------------------------
! Numerical constants                                                         
!-----------------------------------------------------------------------------

  REAL(8), PARAMETER :: ZERO = 0.0
  REAL(8), PARAMETER :: HALF = 0.5
  REAL(8), PARAMETER :: ONE = 1.0
  REAL(8), PARAMETER :: TWO = 2.0


!-----------------------------------------------------------------------------
! Variable species indices                                                    
!-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: IND_SO4 = 0
  INTEGER, PARAMETER :: IND_MSA = 1
  INTEGER, PARAMETER :: IND_CO2 = 2
  INTEGER, PARAMETER :: IND_DRYDEP = 3
  INTEGER, PARAMETER :: IND_N2O = 4
  INTEGER, PARAMETER :: IND_O1D = 5
  INTEGER, PARAMETER :: IND_LISOPOH = 6
  INTEGER, PARAMETER :: IND_CHBr3 = 7
  INTEGER, PARAMETER :: IND_CH2Br2 = 8
  INTEGER, PARAMETER :: IND_CH3Br = 9
  INTEGER, PARAMETER :: IND_H2O2 = 10
  INTEGER, PARAMETER :: IND_PPN = 11
  INTEGER, PARAMETER :: IND_BrNO2 = 12
  INTEGER, PARAMETER :: IND_SO2 = 13
  INTEGER, PARAMETER :: IND_PAN = 14
  INTEGER, PARAMETER :: IND_ALK4 = 15
  INTEGER, PARAMETER :: IND_HNO2 = 16
  INTEGER, PARAMETER :: IND_N2O5 = 17
  INTEGER, PARAMETER :: IND_HNO4 = 18
  INTEGER, PARAMETER :: IND_MAOP = 19
  INTEGER, PARAMETER :: IND_MAP = 20
  INTEGER, PARAMETER :: IND_MP = 21
  INTEGER, PARAMETER :: IND_GLYX = 22
  INTEGER, PARAMETER :: IND_ETP = 23
  INTEGER, PARAMETER :: IND_R4P = 24
  INTEGER, PARAMETER :: IND_RA3P = 25
  INTEGER, PARAMETER :: IND_RB3P = 26
  INTEGER, PARAMETER :: IND_RP = 27
  INTEGER, PARAMETER :: IND_DMS = 28
  INTEGER, PARAMETER :: IND_C3H8 = 29
  INTEGER, PARAMETER :: IND_GP = 30
  INTEGER, PARAMETER :: IND_PP = 31
  INTEGER, PARAMETER :: IND_PRPN = 32
  INTEGER, PARAMETER :: IND_INPN = 33
  INTEGER, PARAMETER :: IND_HOBr = 34
  INTEGER, PARAMETER :: IND_Br2 = 35
  INTEGER, PARAMETER :: IND_HBr = 36
  INTEGER, PARAMETER :: IND_MRP = 37
  INTEGER, PARAMETER :: IND_BrNO3 = 38
  INTEGER, PARAMETER :: IND_C2H6 = 39
  INTEGER, PARAMETER :: IND_ISNP = 40
  INTEGER, PARAMETER :: IND_IAP = 41
  INTEGER, PARAMETER :: IND_VRP = 42
  INTEGER, PARAMETER :: IND_PMN = 43
  INTEGER, PARAMETER :: IND_RIP = 44
  INTEGER, PARAMETER :: IND_ISOP = 45
  INTEGER, PARAMETER :: IND_PRPE = 46
  INTEGER, PARAMETER :: IND_CO = 47
  INTEGER, PARAMETER :: IND_GLYC = 48
  INTEGER, PARAMETER :: IND_ROH = 49
  INTEGER, PARAMETER :: IND_ACET = 50
  INTEGER, PARAMETER :: IND_BrO = 51
  INTEGER, PARAMETER :: IND_MAN2 = 52
  INTEGER, PARAMETER :: IND_HNO3 = 53
  INTEGER, PARAMETER :: IND_A3O2 = 54
  INTEGER, PARAMETER :: IND_R4N1 = 55
  INTEGER, PARAMETER :: IND_MRO2 = 56
  INTEGER, PARAMETER :: IND_RIO1 = 57
  INTEGER, PARAMETER :: IND_IALD = 58
  INTEGER, PARAMETER :: IND_HAC = 59
  INTEGER, PARAMETER :: IND_VRO2 = 60
  INTEGER, PARAMETER :: IND_PRN1 = 61
  INTEGER, PARAMETER :: IND_INO2 = 62
  INTEGER, PARAMETER :: IND_ISN1 = 63
  INTEGER, PARAMETER :: IND_PO2 = 64
  INTEGER, PARAMETER :: IND_B3O2 = 65
  INTEGER, PARAMETER :: IND_ATO2 = 66
  INTEGER, PARAMETER :: IND_IAO2 = 67
  INTEGER, PARAMETER :: IND_GCO3 = 68
  INTEGER, PARAMETER :: IND_KO2 = 69
  INTEGER, PARAMETER :: IND_RCHO = 70
  INTEGER, PARAMETER :: IND_CH2O = 71
  INTEGER, PARAMETER :: IND_MGLY = 72
  INTEGER, PARAMETER :: IND_R4O2 = 73
  INTEGER, PARAMETER :: IND_R4N2 = 74
  INTEGER, PARAMETER :: IND_MAO3 = 75
  INTEGER, PARAMETER :: IND_MVK = 76
  INTEGER, PARAMETER :: IND_RIO2 = 77
  INTEGER, PARAMETER :: IND_MACR = 78
  INTEGER, PARAMETER :: IND_ALD2 = 79
  INTEGER, PARAMETER :: IND_MEK = 80
  INTEGER, PARAMETER :: IND_ETO2 = 81
  INTEGER, PARAMETER :: IND_NO = 82
  INTEGER, PARAMETER :: IND_HO2 = 83
  INTEGER, PARAMETER :: IND_NO3 = 84
  INTEGER, PARAMETER :: IND_Br = 85
  INTEGER, PARAMETER :: IND_NO2 = 86
  INTEGER, PARAMETER :: IND_MCO3 = 87
  INTEGER, PARAMETER :: IND_OH = 88
  INTEGER, PARAMETER :: IND_O3 = 89
  INTEGER, PARAMETER :: IND_RCO3 = 90
  INTEGER, PARAMETER :: IND_MO2 = 91


!-----------------------------------------------------------------------------
! Fixed species indices                                                       
!-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: IND_ACTA = 92
  INTEGER, PARAMETER :: IND_CH4 = 93
  INTEGER, PARAMETER :: IND_EMISSION = 94
  INTEGER, PARAMETER :: IND_EOH = 95
  INTEGER, PARAMETER :: IND_H = 96
  INTEGER, PARAMETER :: IND_H2 = 97
  INTEGER, PARAMETER :: IND_H2O = 98
  INTEGER, PARAMETER :: IND_HCOOH = 99
  INTEGER, PARAMETER :: IND_MOH = 100
  INTEGER, PARAMETER :: IND_N2 = 101
  INTEGER, PARAMETER :: IND_NH2 = 102
  INTEGER, PARAMETER :: IND_NH3 = 103
  INTEGER, PARAMETER :: IND_O = 104
  INTEGER, PARAMETER :: IND_O2 = 105
  INTEGER, PARAMETER :: IND_RCOOH = 106


END MODULE gckpp_parameters
!------------------------ END gckpp_parameters.f90 END -----------------------
