Sat Sep 25 12:40:19 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/TD1/dat7.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m7.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1639.60655084151        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.5959E+01  3.4922E+00 -3.6410E+01  5.8622E+01  6.0528E+01 -3.9179E+01  1.0419E+01  6.8876E+00  1.2687E+01 -6.0734E+00
            -3.1205E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1643.77004301084        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0271E+00  1.0227E+00  1.0527E+00  9.9387E-01  1.0074E+00  1.1086E+00  8.8842E-01  9.3810E-01  9.3404E-01  1.0108E+00
             1.0979E+00
 PARAMETER:  1.2670E-01  1.2248E-01  1.5135E-01  9.3848E-02  1.0734E-01  2.0312E-01 -1.8316E-02  3.6103E-02  3.1765E-02  1.1073E-01
             1.9344E-01
 GRADIENT:   1.2793E+02  2.9476E+01 -1.2797E+01  5.9892E+01  2.6272E+01  7.2593E+00  2.5207E+00  2.6373E+00 -7.5521E+00 -5.9509E+00
             6.1830E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1645.94731244782        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0160E+00  1.1995E+00  9.0080E-01  8.7275E-01  1.0124E+00  1.1161E+00  5.0977E-01  8.2658E-01  1.1942E+00  1.0558E+00
             1.0537E+00
 PARAMETER:  1.1589E-01  2.8192E-01 -4.4683E-03 -3.6110E-02  1.1228E-01  2.0987E-01 -5.7379E-01 -9.0459E-02  2.7744E-01  1.5429E-01
             1.5234E-01
 GRADIENT:   1.0573E+02  5.0797E+01 -1.1178E+01  5.8058E+01  5.3039E+00  1.2507E+01 -8.4761E-01  6.2961E-01  6.0947E+00 -6.9121E+00
            -9.9975E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1648.20949730706        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.6997E-01  1.0932E+00  1.0602E+00  9.1418E-01  1.0372E+00  1.0781E+00  4.9354E-01  8.7779E-01  1.1182E+00  1.1401E+00
             1.0657E+00
 PARAMETER:  6.9506E-02  1.8915E-01  1.5844E-01  1.0274E-02  1.3657E-01  1.7524E-01 -6.0616E-01 -3.0345E-02  2.1171E-01  2.3113E-01
             1.6366E-01
 GRADIENT:   1.7034E+01  8.1209E+00 -4.6194E+00  1.1521E+01  3.7027E+00 -8.5529E-02  1.5064E+00  4.7939E-01  2.3957E+00  1.4471E+00
            -4.1725E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1648.22495205744        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.6563E-01  1.0673E+00  1.1128E+00  9.2758E-01  1.0444E+00  1.0785E+00  4.2454E-01  9.5191E-01  1.1180E+00  1.1501E+00
             1.0697E+00
 PARAMETER:  6.5030E-02  1.6514E-01  2.0692E-01  2.4822E-02  1.4345E-01  1.7558E-01 -7.5674E-01  5.0720E-02  2.1153E-01  2.3984E-01
             1.6734E-01
 GRADIENT:   9.0058E+00  4.0335E+00 -2.8590E+00  6.2868E+00  1.5868E+00  5.4881E-02  1.3648E+00  5.2757E-01  2.4949E+00  1.4330E+00
            -2.4617E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1648.22603704655        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.6353E-01  1.0488E+00  1.1449E+00  9.3762E-01  1.0476E+00  1.0786E+00  3.6904E-01  9.8652E-01  1.1157E+00  1.1550E+00
             1.0722E+00
 PARAMETER:  6.2850E-02  1.4764E-01  2.3536E-01  3.5587E-02  1.4651E-01  1.7570E-01 -8.9685E-01  8.6430E-02  2.0944E-01  2.4414E-01
             1.6969E-01
 GRADIENT:   5.1313E+00  2.0092E+00 -1.8398E+00  3.6421E+00  7.2202E-01  7.9202E-02  1.1157E+00  4.4389E-01  2.1879E+00  1.1911E+00
            -1.5373E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1648.94309437470        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  9.8702E-01  1.1208E+00  1.1713E+00  8.9203E-01  1.0854E+00  1.1164E+00  2.8925E-01  1.1002E+00  1.1772E+00  1.1632E+00
             1.0803E+00
 PARAMETER:  8.6937E-02  2.1408E-01  2.5810E-01 -1.4254E-02  1.8199E-01  2.1011E-01 -1.1405E+00  1.9552E-01  2.6318E-01  2.5120E-01
             1.7724E-01
 GRADIENT:   7.8915E+00  1.0872E+00  1.0071E+00  9.9203E-01 -5.2560E-01  2.4796E+00  1.2317E-01 -1.5587E-01 -1.1671E+00 -8.8355E-01
            -3.3179E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1648.98665797405        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      748
 NPARAMETR:  9.8520E-01  1.1318E+00  1.1617E+00  8.8437E-01  1.0878E+00  1.1120E+00  2.7158E-01  1.1021E+00  1.1908E+00  1.1680E+00
             1.0804E+00
 PARAMETER:  8.5086E-02  2.2383E-01  2.4988E-01 -2.2879E-02  1.8412E-01  2.0616E-01 -1.2035E+00  1.9725E-01  2.7460E-01  2.5527E-01
             1.7730E-01
 GRADIENT:   4.4920E+00  1.2324E+00  1.5764E-01  1.5359E+00  2.6069E-02  9.1160E-01  4.8833E-02 -4.2423E-02 -1.0052E+00 -4.2257E-01
            -2.3966E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1648.99903761738        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      913
 NPARAMETR:  9.8240E-01  1.1320E+00  1.1609E+00  8.8333E-01  1.0878E+00  1.1089E+00  2.6274E-01  1.1052E+00  1.1945E+00  1.1705E+00
             1.0810E+00
 PARAMETER:  8.2239E-02  2.2399E-01  2.4920E-01 -2.4060E-02  1.8412E-01  2.0336E-01 -1.2366E+00  1.9999E-01  2.7771E-01  2.5741E-01
             1.7789E-01
 GRADIENT:  -5.8595E-01 -5.5844E-02 -2.0519E-01  5.7091E-01  4.8585E-02 -1.8081E-01  6.0761E-02  7.8107E-02 -5.9478E-01 -1.4346E-02
             2.0636E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1649.00916924237        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1091
 NPARAMETR:  9.8268E-01  1.1320E+00  1.1691E+00  8.8165E-01  1.0900E+00  1.1094E+00  2.0239E-01  1.1219E+00  1.2085E+00  1.1744E+00
             1.0806E+00
 PARAMETER:  8.2530E-02  2.2399E-01  2.5624E-01 -2.5961E-02  1.8619E-01  2.0384E-01 -1.4975E+00  2.1502E-01  2.8940E-01  2.6073E-01
             1.7754E-01
 GRADIENT:  -1.1118E-02 -6.5006E-02  3.8394E-02 -4.5440E-02 -4.8768E-01 -3.4082E-02 -8.5012E-03 -2.1443E-02 -3.1558E-02 -4.8463E-02
            -2.5059E-02

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1649.00961350472        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1218
 NPARAMETR:  9.8267E-01  1.1320E+00  1.1727E+00  8.8184E-01  1.0917E+00  1.1095E+00  2.0437E-01  1.1278E+00  1.2082E+00  1.1753E+00
             1.0807E+00
 PARAMETER:  8.2522E-02  2.2399E-01  2.5935E-01 -2.5746E-02  1.8774E-01  2.0391E-01 -1.4878E+00  2.2024E-01  2.8913E-01  2.6149E-01
             1.7758E-01
 GRADIENT:  -2.7973E-03 -2.5175E-01 -1.3140E-03 -6.3270E-04 -4.3840E-03 -6.8099E-04  2.4833E-04  1.5214E-03  4.0146E-03 -3.6782E-03
            -5.1926E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1218
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1047E-04 -1.8166E-02 -2.5927E-02 -3.8604E-04 -2.9594E-02
 SE:             2.9873E-02  5.5288E-03  1.2324E-02  2.8549E-02  2.3806E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8103E-01  1.0174E-03  3.5398E-02  9.8921E-01  2.1382E-01

 ETASHRINKSD(%)  1.0000E-10  8.1478E+01  5.8713E+01  4.3559E+00  2.0247E+01
 ETASHRINKVR(%)  1.0000E-10  9.6569E+01  8.2954E+01  8.5221E+00  3.6395E+01
 EBVSHRINKSD(%)  4.0053E-01  8.3735E+01  6.3048E+01  4.4343E+00  1.7329E+01
 EBVSHRINKVR(%)  7.9945E-01  9.7354E+01  8.6346E+01  8.6720E+00  3.1654E+01
 RELATIVEINF(%)  9.9041E+01  2.4631E-01  4.2702E+00  1.0176E+01  2.1064E+01
 EPSSHRINKSD(%)  4.3592E+01
 EPSSHRINKVR(%)  6.8181E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1649.0096135047177     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -913.85878694097948     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.04
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1649.010       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.13E+00  1.17E+00  8.82E-01  1.09E+00  1.11E+00  2.04E-01  1.13E+00  1.21E+00  1.18E+00  1.08E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.23E+02
 
 TH 2
+       -1.37E+01  5.98E+02
 
 TH 3
+        5.44E+00  8.22E+01  9.47E+01
 
 TH 4
+       -1.60E+01  5.99E+02 -1.31E+01  8.94E+02
 
 TH 5
+        1.94E+00 -2.10E+02 -1.52E+02 -2.93E+01  4.68E+02
 
 TH 6
+       -5.46E-01 -3.45E+00  1.67E+00 -1.56E+00 -2.56E-01  1.59E+02
 
 TH 7
+        1.25E+00 -3.57E+01  9.30E-01 -1.67E+01  6.75E-02  5.35E-01  7.30E+00
 
 TH 8
+       -1.47E+00 -1.56E+01 -1.66E+01 -3.22E+00  3.84E+00 -4.23E-01  1.98E+00  8.35E+00
 
 TH 9
+        2.66E+00 -8.98E+01  1.84E+00  1.17E+01  1.38E+00 -1.49E-01  2.06E+01 -2.62E-01  1.12E+02
 
 TH10
+        5.13E-01 -5.26E+00 -8.67E+00 -1.93E+00 -4.61E+01  3.23E-01  4.29E+00  1.03E+01  8.22E-01  6.89E+01
 
 TH11
+       -6.06E+00 -3.61E+01 -1.83E+01 -1.42E+01 -7.45E-01  1.16E+00  3.22E+00  1.02E+01  8.98E+00  1.25E+01  1.82E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       20.066
Stop Time:
Sat Sep 25 12:40:40 CDT 2021
