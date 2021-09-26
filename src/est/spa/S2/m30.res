Sat Sep 25 12:13:14 CDT 2021
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
$DATA ../../../../data/spa/S2/dat30.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.17148111646        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3139E+01 -1.8620E+01 -3.9127E+01  3.6328E+01  3.0797E+01  2.4064E+00 -1.3837E+01  6.3729E+00  5.1592E+00 -8.7494E+00
             9.5463E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1697.73164282038        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0122E+00  1.1761E+00  1.1722E+00  9.0026E-01  1.1530E+00  9.9454E-01  1.2139E+00  9.4609E-01  9.4051E-01  1.1754E+00
             1.0014E+00
 PARAMETER:  1.1210E-01  2.6222E-01  2.5884E-01 -5.0730E-03  2.4238E-01  9.4524E-02  2.9380E-01  4.4586E-02  3.8667E-02  2.6164E-01
             1.0145E-01
 GRADIENT:   5.7449E+01  2.7991E+01 -3.3598E-01  2.5806E+01  1.0002E+01  4.6078E-01  8.0689E+00 -1.8937E+00  2.1294E+00  2.6661E-01
             7.0200E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1698.19010443371        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      265
 NPARAMETR:  1.0214E+00  1.3027E+00  1.0800E+00  8.1672E-01  1.1769E+00  1.0117E+00  1.1080E+00  1.0274E+00  9.3963E-01  1.1783E+00
             9.8547E-01
 PARAMETER:  1.2119E-01  3.6441E-01  1.7696E-01 -1.0246E-01  2.6284E-01  1.1166E-01  2.0259E-01  1.2704E-01  3.7733E-02  2.6408E-01
             8.5368E-02
 GRADIENT:   1.9487E+01  1.5162E+01  2.2774E-01  1.4939E+01  5.4893E+00  4.2547E-02  4.3368E-01 -1.2808E+00 -5.3312E+00 -5.7058E-01
             3.9023E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.73989322150        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.0116E+00  1.3612E+00  1.0269E+00  7.6866E-01  1.1782E+00  1.0115E+00  1.0359E+00  1.1009E+00  1.0317E+00  1.1603E+00
             9.8089E-01
 PARAMETER:  1.1151E-01  4.0840E-01  1.2658E-01 -1.6311E-01  2.6396E-01  1.1145E-01  1.3526E-01  1.9610E-01  1.3121E-01  2.4868E-01
             8.0702E-02
 GRADIENT:  -1.9709E+00  2.3314E+00  4.6183E-01  3.5012E+00 -1.3670E+00  1.8818E-02 -8.1695E-01  6.9935E-02 -6.3179E-01 -1.1232E+00
            -1.5465E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1699.25334524430        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  1.0182E+00  1.6799E+00  6.6058E-01  5.5155E-01  1.2193E+00  1.0130E+00  9.0281E-01  6.7146E-01  1.2612E+00  1.1833E+00
             9.8012E-01
 PARAMETER:  1.1800E-01  6.1873E-01 -3.1463E-01 -4.9502E-01  2.9830E-01  1.1295E-01 -2.2469E-03 -2.9830E-01  3.3210E-01  2.6831E-01
             7.9919E-02
 GRADIENT:   7.9392E+00 -2.8306E+00 -2.4175E+00  2.4775E+00  2.7863E+00 -3.2644E-01  3.3831E-01  3.6423E-01  2.9204E-01  1.4599E+00
             6.1202E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1699.41688211050        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  1.0155E+00  1.8234E+00  5.4142E-01  4.5585E-01  1.2546E+00  1.0144E+00  8.5374E-01  4.2679E-01  1.4218E+00  1.1923E+00
             9.8076E-01
 PARAMETER:  1.1535E-01  7.0072E-01 -5.1356E-01 -6.8560E-01  3.2679E-01  1.1428E-01 -5.8128E-02 -7.5145E-01  4.5194E-01  2.7588E-01
             8.0574E-02
 GRADIENT:   1.3970E+00  6.6011E-01 -2.6335E-01  3.3126E-01 -2.4587E-01  1.0589E-02  2.0403E-01  1.5832E-01  2.5520E-01  2.2731E-01
             3.4737E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1699.47822616181        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  1.0149E+00  1.8150E+00  5.3898E-01  4.6148E-01  1.2485E+00  1.0145E+00  8.5434E-01  1.8837E-01  1.4146E+00  1.1877E+00
             9.8003E-01
 PARAMETER:  1.1479E-01  6.9610E-01 -5.1809E-01 -6.7331E-01  3.2191E-01  1.1438E-01 -5.7421E-02 -1.5694E+00  4.4682E-01  2.7206E-01
             7.9824E-02
 GRADIENT:   2.4882E-01  1.2701E+00 -4.1357E-01  1.2156E+00  9.3623E-01  2.7193E-03 -4.8432E-01  3.0860E-02  1.3456E-01 -8.9073E-02
            -1.6436E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1699.49475548179        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  1.0148E+00  1.8151E+00  5.3706E-01  4.6040E-01  1.2471E+00  1.0144E+00  8.5655E-01  3.2895E-02  1.4145E+00  1.1878E+00
             9.8003E-01
 PARAMETER:  1.1470E-01  6.9612E-01 -5.2164E-01 -6.7566E-01  3.2080E-01  1.1433E-01 -5.4839E-02 -3.3144E+00  4.4674E-01  2.7209E-01
             7.9824E-02
 GRADIENT:   6.9670E-03  4.0517E-01  3.5005E-02  6.1351E-02 -1.2133E-01  1.0371E-03  8.0355E-02  9.6930E-04 -1.6284E-02 -4.2235E-03
             2.0130E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1699.49528283002        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1325
 NPARAMETR:  1.0148E+00  1.8142E+00  5.3750E-01  4.6079E-01  1.2470E+00  1.0144E+00  8.5645E-01  1.1735E-02  1.4140E+00  1.1877E+00
             9.8001E-01
 PARAMETER:  1.1471E-01  6.9565E-01 -5.2083E-01 -6.7481E-01  3.2072E-01  1.1431E-01 -5.4960E-02 -4.3452E+00  4.4642E-01  2.7203E-01
             7.9808E-02
 GRADIENT:   7.0791E-03  3.0111E-02  4.6315E-03  1.6412E-03  7.5256E-03 -2.0420E-03  1.1127E-03  2.8722E-04 -1.0625E-02 -8.5232E-03
            -4.9866E-03

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1699.49528356068        NO. OF FUNC. EVALS.:  67
 CUMULATIVE NO. OF FUNC. EVALS.:     1392
 NPARAMETR:  1.0148E+00  1.8136E+00  5.3737E-01  4.6059E-01  1.2470E+00  1.0145E+00  8.5638E-01  1.0000E-02  1.4146E+00  1.1880E+00
             9.8005E-01
 PARAMETER:  1.1473E-01  6.9575E-01 -5.2099E-01 -6.7500E-01  3.2077E-01  1.1429E-01 -5.5009E-02 -4.5099E+00  4.4654E-01  2.7205E-01
             7.9813E-02
 GRADIENT:   1.0036E-02  2.4047E-01  2.6388E-03  2.8612E-02  1.9442E-02 -1.2427E-02  2.4623E-03  0.0000E+00 -8.6814E-03 -1.0664E-02
            -4.7521E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1392
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.7851E-04 -2.5715E-02 -2.5643E-04  2.7437E-02 -3.6842E-02
 SE:             2.9850E-02  2.5077E-02  8.5127E-05  2.0369E-02  2.3164E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9256E-01  3.0517E-01  2.5931E-03  1.7799E-01  1.1172E-01

 ETASHRINKSD(%)  1.0000E-10  1.5988E+01  9.9715E+01  3.1761E+01  2.2399E+01
 ETASHRINKVR(%)  1.0000E-10  2.9419E+01  9.9999E+01  5.3434E+01  3.9781E+01
 EBVSHRINKSD(%)  3.9995E-01  1.5080E+01  9.9771E+01  3.6294E+01  1.8971E+01
 EBVSHRINKVR(%)  7.9830E-01  2.7887E+01  9.9999E+01  5.9415E+01  3.4343E+01
 RELATIVEINF(%)  9.9074E+01  4.6465E+00  5.8130E-05  2.1698E+00  1.8390E+01
 EPSSHRINKSD(%)  4.4354E+01
 EPSSHRINKVR(%)  6.9035E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1699.4952835606832     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -964.34445699694504     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.31
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1699.495       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.81E+00  5.37E-01  4.61E-01  1.25E+00  1.01E+00  8.56E-01  1.00E-02  1.41E+00  1.19E+00  9.80E-01
 


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
+        1.04E+03
 
 TH 2
+       -5.20E+00  3.47E+02
 
 TH 3
+        9.48E+00  8.15E+01  2.51E+02
 
 TH 4
+       -1.40E+01  3.58E+02 -2.96E+02  1.11E+03
 
 TH 5
+       -3.29E+00 -8.84E+01 -1.84E+02  2.37E+02  3.38E+02
 
 TH 6
+        7.11E+00 -1.06E+00  1.22E+00 -6.19E+00 -8.31E-01  1.97E+02
 
 TH 7
+       -2.71E+00  1.19E+01 -4.72E+00 -2.40E+01 -4.42E+00  3.10E+00  1.67E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.76E+01
 
 TH 9
+        1.59E+00 -1.28E+01 -3.04E+01  5.77E+01  3.21E+00 -1.52E+00  1.48E+01  0.00E+00  2.69E+01
 
 TH10
+        4.41E-03 -1.20E+01 -2.40E+01  1.28E+01 -4.47E+01  1.73E+00  3.20E+00  0.00E+00 -1.69E+01  6.85E+01
 
 TH11
+       -7.33E+00 -1.82E+01 -3.13E+01  2.90E+01  1.02E+00  1.11E+00  1.96E+01  0.00E+00  4.56E+00  1.51E+01  2.27E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.004
Stop Time:
Sat Sep 25 12:13:40 CDT 2021
