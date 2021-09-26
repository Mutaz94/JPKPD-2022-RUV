Sat Sep 25 11:41:55 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1618.06560095993        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9528E+01 -2.2364E+01 -4.2422E+01  3.6747E+01  9.4222E+01  5.4250E+01 -8.8264E+00  1.0537E+00 -9.6992E+00 -1.1179E+01
            -1.1856E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1637.63643743734        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9299E-01  1.0299E+00  1.0443E+00  9.6437E-01  9.5453E-01  8.4662E-01  1.0282E+00  9.8559E-01  1.0323E+00  9.8141E-01
             1.2558E+00
 PARAMETER:  9.2962E-02  1.2945E-01  1.4333E-01  6.3716E-02  5.3460E-02 -6.6503E-02  1.2779E-01  8.5482E-02  1.3175E-01  8.1232E-02
             3.2773E-01
 GRADIENT:  -1.7683E+01 -3.4390E+00  6.0323E+00 -1.1460E+01 -1.4204E+01 -3.8342E+00 -1.4383E+00  1.1603E+00 -1.1545E+00  1.8767E+00
             1.5578E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1638.30952748652        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.9948E-01  1.0962E+00  8.4754E-01  9.2394E-01  9.0404E-01  8.5283E-01  1.1280E+00  7.3511E-01  1.0303E+00  8.8374E-01
             1.2539E+00
 PARAMETER:  9.9475E-02  1.9184E-01 -6.5412E-02  2.0890E-02 -8.7996E-04 -5.9195E-02  2.2047E-01 -2.0773E-01  1.2986E-01 -2.3597E-02
             3.2628E-01
 GRADIENT:  -2.4463E+00  8.4241E+00  1.3506E+00  5.7910E+00 -5.1847E+00 -1.2462E+00  3.6800E+00  8.3912E-01  3.4126E+00 -2.0974E+00
             1.6380E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1638.57864935354        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      250
 NPARAMETR:  1.0004E+00  1.1120E+00  8.3359E-01  9.1203E-01  9.1178E-01  8.5504E-01  1.0767E+00  5.9404E-01  1.0388E+00  9.3590E-01
             1.2557E+00
 PARAMETER:  1.0043E-01  2.0614E-01 -8.2014E-02  7.9173E-03  7.6397E-03 -5.6605E-02  1.7395E-01 -4.2080E-01  1.3804E-01  3.3758E-02
             3.2768E-01
 GRADIENT:  -2.8597E+01 -2.7172E+00 -1.7358E+00  9.6372E-01 -1.7597E+00 -2.7170E+00 -1.7926E-01  4.1238E-01  9.9647E-01  1.7051E+00
             2.1523E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1639.22491093555        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      428
 NPARAMETR:  1.0137E+00  1.3848E+00  6.7952E-01  7.3646E-01  9.7433E-01  8.6481E-01  9.1752E-01  3.3970E-01  1.2099E+00  9.5506E-01
             1.2531E+00
 PARAMETER:  1.1358E-01  4.2552E-01 -2.8637E-01 -2.0590E-01  7.3994E-02 -4.5241E-02  1.3915E-02 -9.7969E-01  2.9052E-01  5.4014E-02
             3.2564E-01
 GRADIENT:   7.2970E+00  4.9287E-01 -1.0420E+00  3.2684E+00  2.1214E-01  1.4925E+00  7.9448E-02  8.0635E-02  8.2118E-01  6.3810E-01
             2.5161E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1639.36669244614        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  1.0114E+00  1.5447E+00  5.9655E-01  6.2950E-01  1.0227E+00  8.6170E-01  8.4792E-01  1.8266E-01  1.3444E+00  9.7317E-01
             1.2544E+00
 PARAMETER:  1.1138E-01  5.3480E-01 -4.1659E-01 -3.6283E-01  1.2242E-01 -4.8844E-02 -6.4965E-02 -1.6001E+00  3.9597E-01  7.2801E-02
             3.2669E-01
 GRADIENT:   3.0351E-01  2.2960E-01  1.2379E-02  2.8865E-03 -1.1076E-01  8.0601E-02  9.0802E-02  2.5722E-02 -1.4440E-02  9.0581E-02
             2.3031E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1639.36964144761        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  1.0112E+00  1.5566E+00  5.8577E-01  6.2114E-01  1.0251E+00  8.6125E-01  8.4298E-01  1.3666E-01  1.3555E+00  9.7239E-01
             1.2544E+00
 PARAMETER:  1.1111E-01  5.4248E-01 -4.3483E-01 -3.7621E-01  1.2474E-01 -4.9367E-02 -7.0816E-02 -1.8902E+00  4.0419E-01  7.2002E-02
             3.2665E-01
 GRADIENT:  -6.0785E-01 -9.5676E-01 -6.4350E-01  2.1669E-01  1.2369E+00 -1.3291E-01 -4.8190E-02  1.6301E-02  5.4653E-02  4.6868E-02
             9.1524E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1639.37877795367        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      957
 NPARAMETR:  1.0113E+00  1.5441E+00  5.9302E-01  6.2944E-01  1.0206E+00  8.6150E-01  8.4846E-01  2.8040E-02  1.3436E+00  9.7160E-01
             1.2545E+00
 PARAMETER:  1.1126E-01  5.3445E-01 -4.2254E-01 -3.6293E-01  1.2043E-01 -4.9082E-02 -6.4337E-02 -3.4741E+00  3.9533E-01  7.1187E-02
             3.2674E-01
 GRADIENT:  -8.6594E-02 -1.5525E-02  4.2741E-02 -7.4013E-02 -1.1045E-01 -1.6841E-02  7.0149E-03  5.5443E-04  1.8836E-02  1.5707E-02
             1.9329E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1639.37904724647        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1134
 NPARAMETR:  1.0114E+00  1.5433E+00  5.9301E-01  6.2995E-01  1.0203E+00  8.6154E-01  8.4883E-01  1.0000E-02  1.3424E+00  9.7118E-01
             1.2544E+00
 PARAMETER:  1.1130E-01  5.3395E-01 -4.2255E-01 -3.6211E-01  1.2006E-01 -4.9034E-02 -6.3895E-02 -4.9258E+00  3.9449E-01  7.0753E-02
             3.2667E-01
 GRADIENT:   1.4760E-02 -4.3219E-02 -1.4572E-02 -1.6506E-03  3.7559E-02 -4.8282E-04 -6.2623E-03  0.0000E+00 -1.3020E-03 -5.2733E-03
            -3.7396E-03

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1639.37904724647        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1156
 NPARAMETR:  1.0114E+00  1.5433E+00  5.9301E-01  6.2995E-01  1.0203E+00  8.6154E-01  8.4883E-01  1.0000E-02  1.3424E+00  9.7118E-01
             1.2544E+00
 PARAMETER:  1.1130E-01  5.3395E-01 -4.2255E-01 -3.6211E-01  1.2006E-01 -4.9034E-02 -6.3895E-02 -4.9258E+00  3.9449E-01  7.0753E-02
             3.2667E-01
 GRADIENT:   1.4760E-02 -4.3219E-02 -1.4572E-02 -1.6506E-03  3.7559E-02 -4.8282E-04 -6.2623E-03  0.0000E+00 -1.3020E-03 -5.2733E-03
            -3.7396E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1156
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.6861E-05 -2.5217E-02 -2.5609E-04  1.7861E-02 -3.1394E-02
 SE:             2.9698E-02  2.2504E-02  1.0419E-04  2.2940E-02  2.2267E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9874E-01  2.6247E-01  1.3972E-02  4.3623E-01  1.5857E-01

 ETASHRINKSD(%)  5.0793E-01  2.4608E+01  9.9651E+01  2.3147E+01  2.5403E+01
 ETASHRINKVR(%)  1.0133E+00  4.3161E+01  9.9999E+01  4.0937E+01  4.4353E+01
 EBVSHRINKSD(%)  8.7596E-01  2.3947E+01  9.9695E+01  2.4545E+01  2.4384E+01
 EBVSHRINKVR(%)  1.7443E+00  4.2159E+01  9.9999E+01  4.3065E+01  4.2822E+01
 RELATIVEINF(%)  9.8088E+01  3.0982E+00  9.2331E-05  3.2458E+00  9.6857E+00
 EPSSHRINKSD(%)  4.2448E+01
 EPSSHRINKVR(%)  6.6877E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1639.3790472464716     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -904.22822068273342     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1639.379       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.54E+00  5.93E-01  6.30E-01  1.02E+00  8.62E-01  8.49E-01  1.00E-02  1.34E+00  9.71E-01  1.25E+00
 


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
+        1.44E+03
 
 TH 2
+       -1.16E+01  3.93E+02
 
 TH 3
+        1.44E+01  1.45E+02  3.12E+02
 
 TH 4
+       -2.79E+01  3.36E+02 -2.38E+02  9.00E+02
 
 TH 5
+       -5.78E+00 -2.26E+02 -3.54E+02  2.48E+02  6.45E+02
 
 TH 6
+       -3.29E+00 -2.55E+00  6.95E-01 -4.91E+00 -1.96E+00  2.60E+02
 
 TH 7
+        3.60E+00  1.15E+01  3.37E+00 -1.54E+01 -1.45E+01 -8.21E-01  9.40E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.62E+00 -1.75E+01 -2.56E+01  4.68E+01  2.68E+00  3.80E-01  2.01E+01  0.00E+00  4.30E+01
 
 TH10
+       -1.96E+00 -1.21E+01 -3.16E+01 -6.13E+00 -5.56E+01  3.58E+00  1.42E+01  0.00E+00  5.01E+00  7.94E+01
 
 TH11
+       -1.10E+01 -1.76E+01 -2.86E+01  1.36E+00 -1.93E+00  2.88E+00  1.01E+01  0.00E+00  5.63E+00  1.67E+01  1.39E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.288
Stop Time:
Sat Sep 25 11:42:16 CDT 2021
