Sat Sep 25 01:21:58 CDT 2021
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
$DATA ../../../../data/int/SL2/dat56.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      998
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      898
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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -689.400012374666        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9862E+01 -1.1581E+02  2.9128E+02  6.9595E+01  1.5965E+02 -1.8514E+01 -6.7090E+01 -3.5254E+02 -1.2227E+02 -6.7782E+00
            -5.8465E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2802.32536498030        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0350E+00  1.4419E+00  7.8440E-01  8.1101E-01  1.1127E+00  1.0705E+00  1.0223E+00  1.0291E+00  1.0137E+00  8.3024E-01
             2.5844E+00
 PARAMETER:  1.3437E-01  4.6595E-01 -1.4283E-01 -1.0947E-01  2.0678E-01  1.6808E-01  1.2205E-01  1.2867E-01  1.1365E-01 -8.6045E-02
             1.0495E+00
 GRADIENT:   2.6445E+01  3.1873E+01  8.3645E-01  1.8898E+01 -3.3525E+01  1.0448E+01  6.5403E+00  2.9617E+00 -9.7187E+00 -8.2133E+00
            -3.9460E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2805.52644129011        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0506E+00  1.6561E+00  7.6955E-01  7.2745E-01  1.2377E+00  1.0468E+00  8.6438E-01  7.5237E-01  1.2208E+00  1.0366E+00
             2.5938E+00
 PARAMETER:  1.4941E-01  6.0445E-01 -1.6195E-01 -2.1821E-01  3.1329E-01  1.4573E-01 -4.5743E-02 -1.8453E-01  2.9949E-01  1.3594E-01
             1.0531E+00
 GRADIENT:   5.6545E+01  9.3587E+01  5.6673E+00  6.4449E+01 -4.1021E+01  1.1777E+00 -4.5031E+00 -8.0741E-01  2.1137E+00  4.7559E+00
            -2.6171E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2816.90538791408        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0215E+00  1.9107E+00  5.9272E-01  4.8206E-01  1.4869E+00  1.0524E+00  8.1729E-01  5.7415E-01  1.3385E+00  1.1378E+00
             2.5887E+00
 PARAMETER:  1.2125E-01  7.4748E-01 -4.2303E-01 -6.2969E-01  4.9672E-01  1.5112E-01 -1.0177E-01 -4.5487E-01  3.9156E-01  2.2912E-01
             1.0511E+00
 GRADIENT:   1.4446E+00 -2.3342E+01  5.4250E-01 -1.8559E+00  4.0320E+00  4.9201E+00 -2.5778E+00 -3.9195E-01 -4.3012E+00 -3.2032E+00
            -6.3241E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2823.37508763022        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0178E+00  2.2432E+00  3.3385E-01  2.6542E-01  1.7104E+00  1.0311E+00  7.5403E-01  4.0128E-01  2.0169E+00  1.2922E+00
             2.5757E+00
 PARAMETER:  1.1765E-01  9.0791E-01 -9.9706E-01 -1.2264E+00  6.3672E-01  1.3061E-01 -1.8233E-01 -8.1308E-01  8.0157E-01  3.5632E-01
             1.0461E+00
 GRADIENT:  -6.0026E+00 -8.5809E+00 -2.1133E+00  3.1778E+00  4.1851E+00 -3.1928E+00  3.5760E-01  7.1414E-02  2.6623E+00  1.5314E+00
             4.2178E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2825.29141588033        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      462
 NPARAMETR:  1.0300E+00  2.4242E+00  2.2102E-01  1.7090E-01  1.8339E+00  1.0524E+00  7.3395E-01  2.6164E-01  2.4950E+00  1.3493E+00
             2.5731E+00
 PARAMETER:  1.2954E-01  9.8551E-01 -1.4095E+00 -1.6667E+00  7.0646E-01  1.5107E-01 -2.0931E-01 -1.2408E+00  1.0143E+00  3.9958E-01
             1.0451E+00
 GRADIENT:   9.2291E+00  9.4458E+00 -1.3903E+00  4.6226E+00  5.9154E+00  3.3103E+00  9.1007E-01  4.2395E-02  7.5858E-01  6.5866E-02
             1.3034E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2825.56207504219        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      638
 NPARAMETR:  1.0250E+00  2.4652E+00  1.9512E-01  1.3662E-01  1.8635E+00  1.0426E+00  7.2553E-01  2.3189E-01  2.7644E+00  1.3671E+00
             2.5716E+00
 PARAMETER:  1.2466E-01  1.0023E+00 -1.5341E+00 -1.8905E+00  7.2243E-01  1.4170E-01 -2.2086E-01 -1.3615E+00  1.1168E+00  4.1271E-01
             1.0445E+00
 GRADIENT:   9.4748E-02  4.6663E-01 -4.3086E-02 -8.4633E-03 -1.4656E-01  8.1202E-02  1.1334E-01  3.1835E-02 -8.4007E-04  1.7285E-01
             6.6558E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2825.57472147426        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  1.0247E+00  2.4606E+00  1.9972E-01  1.3968E-01  1.8602E+00  1.0422E+00  7.2559E-01  1.1293E-01  2.7436E+00  1.3641E+00
             2.5709E+00
 PARAMETER:  1.2442E-01  1.0004E+00 -1.5109E+00 -1.8684E+00  7.2068E-01  1.4138E-01 -2.2077E-01 -2.0810E+00  1.1093E+00  4.1053E-01
             1.0443E+00
 GRADIENT:  -3.9021E-01  4.6829E-01  3.0260E-02  9.0095E-02 -5.0423E-02 -5.6705E-02 -1.0626E-01  7.5130E-03  3.1082E-02  5.5406E-02
             9.4601E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2825.57881964008        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      990
 NPARAMETR:  1.0249E+00  2.4624E+00  1.9752E-01  1.3824E-01  1.8618E+00  1.0424E+00  7.2560E-01  1.0946E-02  2.7540E+00  1.3649E+00
             2.5708E+00
 PARAMETER:  1.2462E-01  1.0011E+00 -1.5219E+00 -1.8787E+00  7.2152E-01  1.4152E-01 -2.2075E-01 -4.4148E+00  1.1130E+00  4.1105E-01
             1.0442E+00
 GRADIENT:   2.5762E-02 -2.6743E-02 -3.6118E-03 -4.4114E-03  3.8005E-03  4.2707E-03  7.8675E-03  7.1025E-05 -1.6262E-03 -3.7128E-03
            -1.2586E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -2825.57883242450        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1088
 NPARAMETR:  1.0250E+00  2.4611E+00  1.9745E-01  1.3837E-01  1.8615E+00  1.0425E+00  7.2575E-01  1.0000E-02  2.7531E+00  1.3648E+00
             2.5707E+00
 PARAMETER:  1.2460E-01  1.0011E+00 -1.5209E+00 -1.8779E+00  7.2145E-01  1.4151E-01 -2.2075E-01 -4.9993E+00  1.1128E+00  4.1099E-01
             1.0442E+00
 GRADIENT:  -9.5248E-03  5.9563E-02  9.0827E-04 -2.5030E-04  1.7442E-03 -1.4897E-03 -2.3216E-03  0.0000E+00  6.1674E-05  1.8977E-04
             3.2836E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1088
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4717E-03 -1.9026E-02 -2.1449E-05  3.0412E-02 -1.9509E-02
 SE:             2.9514E-02  2.7118E-02  1.9978E-05  1.5403E-02  2.5836E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6023E-01  4.8293E-01  2.8299E-01  4.8335E-02  4.5019E-01

 ETASHRINKSD(%)  1.1242E+00  9.1510E+00  9.9933E+01  4.8399E+01  1.3446E+01
 ETASHRINKVR(%)  2.2357E+00  1.7465E+01  1.0000E+02  7.3373E+01  2.5084E+01
 EBVSHRINKSD(%)  1.3864E+00  8.7873E+00  9.9921E+01  5.7366E+01  1.0996E+01
 EBVSHRINKVR(%)  2.7535E+00  1.6802E+01  1.0000E+02  8.1823E+01  2.0782E+01
 RELATIVEINF(%)  9.7197E+01  2.3579E+01  3.0495E-05  3.9171E+00  4.4753E+01
 EPSSHRINKSD(%)  1.6457E+01
 EPSSHRINKVR(%)  3.0205E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          898
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1650.4136056355921     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2825.5788324245023     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1175.1652267889101     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.56
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2825.579       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.46E+00  1.98E-01  1.38E-01  1.86E+00  1.04E+00  7.26E-01  1.00E-02  2.75E+00  1.36E+00  2.57E+00
 


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
+        9.50E+02
 
 TH 2
+       -9.05E+00  2.97E+02
 
 TH 3
+        2.57E+00  3.77E+01  2.24E+02
 
 TH 4
+       -2.43E+01  3.31E+02 -3.55E+02  1.78E+03
 
 TH 5
+       -3.24E+00 -2.62E+01 -3.49E+01  1.50E+02  1.37E+02
 
 TH 6
+        5.59E+00 -2.78E+00  8.94E-01 -9.25E+00 -1.21E+00  1.76E+02
 
 TH 7
+        1.17E+00  1.60E+00 -1.00E+01 -2.84E+01 -8.06E-02 -2.44E+00  2.78E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.01E-01 -3.52E+00 -1.40E+01  6.01E+01 -8.22E-01  1.12E-01  1.16E+00  0.00E+00  5.71E+00
 
 TH10
+       -8.33E-02 -3.34E+00  1.22E+01  2.42E+01 -8.87E+00  7.86E-01 -2.07E+00  0.00E+00  1.23E+00  6.51E+01
 
 TH11
+       -1.27E+01 -1.23E+01  3.64E-01 -8.87E+00  4.94E-01  2.07E+00  1.04E+01  0.00E+00  6.08E-01  6.65E+00  1.77E+02
 
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
 #CPUT: Total CPU Time in Seconds,       40.372
Stop Time:
Sat Sep 25 01:22:40 CDT 2021
