Wed Sep 29 16:38:28 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat51.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 NO. OF SIG. FIGURES REQUIRED:            2
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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1549.44109722550        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1915E+02  1.8392E+01 -2.2930E+01  1.0055E+02  1.7693E+02  6.0452E+01  7.9459E+00 -2.4750E+01  4.3558E+00 -6.4901E+01
            -8.7828E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1569.57543946806        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.6069E-01  9.3676E-01  8.8993E-01  9.8059E-01  8.6209E-01  8.7646E-01  9.4221E-01  1.1018E+00  9.7780E-01  1.0564E+00
             1.1424E+00
 PARAMETER:  5.9897E-02  3.4671E-02 -1.6607E-02  8.0404E-02 -4.8397E-02 -3.1863E-02  4.0478E-02  1.9693E-01  7.7551E-02  1.5487E-01
             2.3316E-01
 GRADIENT:   3.2597E+02 -2.4943E+01 -1.2690E+01 -9.5952E+00  5.8467E+01  5.9407E-01  5.9306E+00 -4.8565E+00 -3.6590E-02 -7.7781E+00
            -1.0133E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1575.40050424643        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      242
 NPARAMETR:  9.6383E-01  8.7541E-01  9.9004E-01  1.0523E+00  8.6001E-01  9.3268E-01  6.8608E-01  1.1439E+00  1.0210E+00  1.1663E+00
             1.1687E+00
 PARAMETER:  6.3160E-02 -3.3062E-02  8.9993E-02  1.5102E-01 -5.0811E-02  3.0309E-02 -2.7676E-01  2.3447E-01  1.2081E-01  2.5384E-01
             2.5593E-01
 GRADIENT:   8.4704E+00 -1.3756E+01 -4.4080E+00 -1.8917E+01  1.4166E+01 -1.6585E+01  3.5309E+00 -7.1557E+00  2.4918E+00  4.6726E+00
             1.4470E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1577.33759005408        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      428
 NPARAMETR:  9.6213E-01  9.9606E-01  1.1844E+00  9.8284E-01  9.7634E-01  9.3185E-01  4.3491E-01  1.5834E+00  1.1719E+00  1.2368E+00
             1.1815E+00
 PARAMETER:  6.1395E-02  9.6050E-02  2.6921E-01  8.2689E-02  7.6054E-02  2.9418E-02 -7.3263E-01  5.5958E-01  2.5860E-01  3.1249E-01
             2.6682E-01
 GRADIENT:   5.7658E+00 -1.7637E+01 -2.0565E+00 -9.9818E+00  1.2991E+01 -1.6823E+01  3.5668E+00 -3.3282E+00  1.4333E+01  2.5335E+00
             4.3581E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1577.71697577869        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      621
 NPARAMETR:  9.6189E-01  9.9755E-01  1.1832E+00  9.8235E-01  9.7643E-01  9.7147E-01  4.3497E-01  1.5845E+00  1.1719E+00  1.2359E+00
             1.1808E+00
 PARAMETER:  6.1148E-02  9.7551E-02  2.6825E-01  8.2191E-02  7.6151E-02  7.1058E-02 -7.3247E-01  5.6030E-01  2.5862E-01  3.1176E-01
             2.6618E-01
 GRADIENT:   4.8189E+00 -1.7273E+01 -2.0042E+00 -9.5186E+00  1.2933E+01  2.4354E-03  3.5724E+00 -3.2771E+00  1.4119E+01  2.3875E+00
             4.2985E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1578.56892101363        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  9.5695E-01  1.0021E+00  1.1832E+00  9.8318E-01  9.7132E-01  9.6751E-01  3.5428E-01  1.5858E+00  1.1517E+00  1.2290E+00
             1.1764E+00
 PARAMETER:  5.5999E-02  1.0211E-01  2.6825E-01  8.3039E-02  7.0900E-02  6.6972E-02 -9.3765E-01  5.6107E-01  2.4126E-01  3.0623E-01
             2.6248E-01
 GRADIENT:  -7.5697E+00 -4.8581E+00  4.3690E-01 -2.1385E+00  7.1469E+00 -1.7468E+00  2.1229E+00 -4.1915E+00  5.2025E+00  1.6942E+00
             2.0577E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1579.47546706691        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      928
 NPARAMETR:  9.5783E-01  9.9934E-01  1.1623E+00  9.8406E-01  9.6613E-01  9.7358E-01  1.0271E-01  1.5813E+00  1.1336E+00  1.2227E+00
             1.1727E+00
 PARAMETER:  5.6910E-02  9.9342E-02  2.5041E-01  8.3930E-02  6.5544E-02  7.3225E-02 -2.1758E+00  5.5828E-01  2.2544E-01  3.0110E-01
             2.5931E-01
 GRADIENT:   3.1080E+02  3.9015E+01  8.0365E-01  4.7694E+01  1.5374E+01  4.5981E+01  8.8360E-01 -2.0676E+00  6.9591E+00  4.1243E+00
             2.2611E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1579.57530491223        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1089
 NPARAMETR:  9.5783E-01  9.9934E-01  1.1682E+00  9.8274E-01  9.6613E-01  9.7228E-01  9.2257E-02  1.5813E+00  1.1502E+00  1.2227E+00
             1.1727E+00
 PARAMETER:  5.6910E-02  9.9342E-02  2.5547E-01  8.2592E-02  6.5544E-02  7.1889E-02 -2.2832E+00  5.5828E-01  2.3995E-01  3.0110E-01
             2.5931E-01
 GRADIENT:  -5.6499E+00  6.1144E-01 -6.7045E-01  7.4400E-01  7.8410E+00 -3.8364E-02  1.5094E-01 -4.4548E+00 -2.7057E+00  8.6194E-01
             5.4011E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1579.72207596120        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1245            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6071E-01  9.9911E-01  1.1654E+00  9.8123E-01  9.5668E-01  9.7294E-01  1.0000E-02  1.5857E+00  1.1609E+00  1.2173E+00
             1.1773E+00
 PARAMETER:  5.9919E-02  9.9113E-02  2.5306E-01  8.1056E-02  5.5710E-02  7.2570E-02 -4.6481E+00  5.6105E-01  2.4917E-01  2.9660E-01
             2.6321E-01
 GRADIENT:   3.1575E+02  3.5134E+01  4.9220E+00  4.3839E+01  5.8179E+00  4.5102E+01  0.0000E+00 -2.2040E+00  1.5920E+01  4.4335E+00
             4.4915E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1579.72207596120        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:     1305
 NPARAMETR:  9.6071E-01  9.9911E-01  1.1654E+00  9.8123E-01  9.5668E-01  9.7294E-01  1.0000E-02  1.5857E+00  1.1609E+00  1.2173E+00
             1.1773E+00
 PARAMETER:  5.9919E-02  9.9113E-02  2.5306E-01  8.1056E-02  5.5710E-02  7.2570E-02 -4.6481E+00  5.6105E-01  2.4917E-01  2.9660E-01
             2.6321E-01
 GRADIENT:   1.3020E+00  4.2641E-01  2.0728E+00 -6.7182E-01 -9.1546E-02  2.7801E-01  0.0000E+00 -4.2081E+00 -9.6778E-02  1.4222E+00
             2.3921E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1305
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2911E-03 -8.1196E-04 -3.8719E-02 -2.3440E-03 -3.3078E-02
 SE:             2.9811E-02  2.1936E-04  1.6217E-02  2.8872E-02  2.2805E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6545E-01  2.1434E-04  1.6962E-02  9.3530E-01  1.4693E-01

 ETASHRINKSD(%)  1.3103E-01  9.9265E+01  4.5670E+01  3.2753E+00  2.3601E+01
 ETASHRINKVR(%)  2.6190E-01  9.9995E+01  7.0483E+01  6.4432E+00  4.1632E+01
 EBVSHRINKSD(%)  6.0845E-01  9.9329E+01  5.0793E+01  3.6412E+00  2.0307E+01
 EBVSHRINKVR(%)  1.2132E+00  9.9995E+01  7.5787E+01  7.1497E+00  3.6490E+01
 RELATIVEINF(%)  9.8552E+01  4.4662E-04  6.6973E+00  1.0879E+01  1.9292E+01
 EPSSHRINKSD(%)  4.4836E+01
 EPSSHRINKVR(%)  6.9569E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1579.7220759611964     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -844.57124939745825     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.92
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1579.722       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  9.99E-01  1.17E+00  9.81E-01  9.57E-01  9.73E-01  1.00E-02  1.59E+00  1.16E+00  1.22E+00  1.18E+00
 


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
+        4.50E+07
 
 TH 2
+       -2.08E+01  4.16E+07
 
 TH 3
+        3.84E+02 -1.41E+07  4.78E+06
 
 TH 4
+       -1.76E+01  4.24E+07 -1.44E+07  4.31E+07
 
 TH 5
+        4.52E+07  2.61E+03 -1.47E+07  2.84E+03  4.54E+07
 
 TH 6
+        1.66E-01  1.40E+03  4.77E+02  1.43E+03 -5.55E-01  2.05E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -4.87E+06 -4.68E+06 -1.59E+06 -4.77E+06  1.52E+00 -8.29E-01  0.00E+00  1.05E+06
 
 TH 9
+        1.50E+07  1.44E+07 -4.87E+06  1.03E+03  1.23E+00  7.17E-01  0.00E+00 -1.62E+06  1.27E+02
 
 TH10
+        3.10E+02  1.15E+07  7.81E+06  1.17E+07  7.27E+02  3.89E+02  0.00E+00 -4.74E+03  2.77E+02  3.19E+06
 
 TH11
+       -1.40E+07 -1.35E+07 -3.93E+04 -1.37E+07 -1.41E+07 -1.38E+07  0.00E+00  1.50E+06 -4.66E+06 -3.21E+04  4.37E+06
 
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
 #CPUT: Total CPU Time in Seconds,       22.404
Stop Time:
Wed Sep 29 16:38:53 CDT 2021
