Wed Sep 29 13:37:57 CDT 2021
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
$DATA ../../../../data/spa/A3/dat55.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -103.622808301138        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2880E+02  4.2931E+01  7.6235E+01  2.2217E+00  2.4479E+02  7.6125E+01 -5.7075E+01 -4.2908E+01 -8.4544E+01 -1.9526E+02
            -2.7116E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1080.46752528140        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0222E+00  9.0000E-01  8.1176E-01  1.1009E+00  8.4015E-01  8.6011E-01  1.0957E+00  9.3285E-01  1.1628E+00  1.1820E+00
             1.8379E+00
 PARAMETER:  1.2198E-01 -5.3556E-03 -1.0855E-01  1.9614E-01 -7.4180E-02 -5.0697E-02  1.9136E-01  3.0494E-02  2.5080E-01  2.6719E-01
             7.0865E-01
 GRADIENT:   2.2864E+02  2.5686E+00 -1.5728E+01  4.9174E+01  1.0613E+02  1.5056E+01 -1.1597E+01  2.5055E+00 -2.4866E+00 -3.7914E+01
            -6.5849E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1115.11361814473        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      186
 NPARAMETR:  1.0292E+00  6.2415E-01  2.5539E-01  1.1995E+00  3.1873E-01  7.8217E-01  1.2011E+00  3.7424E-01  1.0113E+00  4.6594E-01
             1.8133E+00
 PARAMETER:  1.2880E-01 -3.7137E-01 -1.2650E+00  2.8193E-01 -1.0434E+00 -1.4568E-01  2.8323E-01 -8.8286E-01  1.1121E-01 -6.6370E-01
             6.9517E-01
 GRADIENT:   1.1085E+02  1.4122E+02  7.2771E+01  2.1366E+02 -1.0879E+02 -3.6786E+01 -3.0140E+01 -8.9559E+00 -6.1813E+01 -2.7135E+01
            -6.1049E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1244.83543789368        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      363
 NPARAMETR:  1.0000E+00  6.3027E-01  1.7840E-01  9.9401E-01  3.1685E-01  8.9672E-01  5.3040E-01  1.9359E-02  1.3872E+00  6.8253E-01
             2.4499E+00
 PARAMETER:  1.0001E-01 -3.6161E-01 -1.6237E+00  9.3991E-02 -1.0493E+00 -9.0147E-03 -5.3413E-01 -3.8446E+00  4.2729E-01 -2.8195E-01
             9.9603E-01
 GRADIENT:   3.8295E+01 -9.7653E+01 -7.4582E+01  1.9971E+01  2.3465E+02  2.1480E+01 -6.8899E+00 -5.3003E-03  2.5672E+01 -1.3801E+01
            -2.1791E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1295.96691576129        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      541
 NPARAMETR:  1.0260E+00  4.0429E-01  1.8415E-01  1.0576E+00  2.3729E-01  8.1574E-01  1.2434E+00  1.5014E-02  1.1065E+00  4.0440E-01
             3.3505E+00
 PARAMETER:  1.2563E-01 -8.0561E-01 -1.5920E+00  1.5605E-01 -1.3385E+00 -1.0366E-01  3.1781E-01 -4.0988E+00  2.0122E-01 -8.0534E-01
             1.3091E+00
 GRADIENT:   2.4425E+01  1.2794E+01  9.1116E+00  3.0478E+00 -2.9160E+01 -5.8553E+00 -6.4152E-01 -3.5500E-03 -3.5007E-01 -4.2040E+00
             7.6618E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1297.44574385288        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      718
 NPARAMETR:  1.0205E+00  3.0319E-01  2.0528E-01  1.1073E+00  2.2963E-01  8.2502E-01  1.6765E+00  1.0000E-02  1.0498E+00  5.4379E-01
             3.2905E+00
 PARAMETER:  1.2034E-01 -1.0934E+00 -1.4834E+00  2.0193E-01 -1.3713E+00 -9.2350E-02  6.1673E-01 -5.2895E+00  1.4862E-01 -5.0920E-01
             1.2910E+00
 GRADIENT:  -4.0706E+00  9.7341E+00  9.9081E+00 -7.3870E-01 -2.7584E+01  8.7354E-01  1.7836E+00  0.0000E+00 -2.8288E+00  2.1011E+00
             5.7415E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1298.46524379273        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      894
 NPARAMETR:  1.0185E+00  2.2796E-01  2.3526E-01  1.1700E+00  2.3901E-01  8.1488E-01  2.1306E+00  1.0000E-02  1.0148E+00  5.5142E-01
             3.3023E+00
 PARAMETER:  1.1833E-01 -1.3786E+00 -1.3471E+00  2.5698E-01 -1.3313E+00 -1.0472E-01  8.5642E-01 -6.5134E+00  1.1471E-01 -4.9526E-01
             1.2946E+00
 GRADIENT:  -2.8063E+00  1.3168E+00  8.9722E+00 -3.6183E+00 -1.4615E+01  5.4517E-01 -1.5603E+00  0.0000E+00  1.1939E+00 -6.9287E-01
            -1.0514E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1298.86441436864        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  1.0158E+00  1.7533E-01  2.6455E-01  1.2302E+00  2.5342E-01  8.0246E-01  2.7478E+00  1.0000E-02  9.6723E-01  5.6761E-01
             3.3351E+00
 PARAMETER:  1.1570E-01 -1.6411E+00 -1.2297E+00  3.0714E-01 -1.2727E+00 -1.2007E-01  1.1108E+00 -7.9346E+00  6.6685E-02 -4.6631E-01
             1.3045E+00
 GRADIENT:   8.4377E-01 -9.0232E-02  1.6308E+00  1.6948E+00 -1.7382E+00 -1.9084E+00 -6.1146E-01  0.0000E+00  6.2563E-02  2.6182E-01
            -3.1933E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1298.88590994110        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  1.0155E+00  1.7106E-01  2.5936E-01  1.2228E+00  2.4987E-01  8.0757E-01  2.7758E+00  1.0000E-02  9.7157E-01  5.6544E-01
             3.3352E+00
 PARAMETER:  1.1531E-01 -1.6640E+00 -1.2496E+00  3.0145E-01 -1.2866E+00 -1.1410E-01  1.1210E+00 -8.0674E+00  7.2154E-02 -4.7092E-01
             1.3038E+00
 GRADIENT:  -2.1971E-01  1.9007E-01 -8.5947E-02  3.5797E-01  3.5015E-01 -8.0654E-02  8.4900E-03  0.0000E+00  7.1022E-02 -3.6994E-02
            -3.7043E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1202
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6319E-04  3.5244E-02 -8.1223E-05 -2.1440E-02  1.3748E-02
 SE:             2.8125E-02  1.4470E-02  2.1036E-04  2.4787E-02  1.6843E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8402E-01  1.4868E-02  6.9941E-01  3.8706E-01  4.1436E-01

 ETASHRINKSD(%)  5.7775E+00  5.1522E+01  9.9295E+01  1.6960E+01  4.3574E+01
 ETASHRINKVR(%)  1.1221E+01  7.6499E+01  9.9995E+01  3.1043E+01  6.8161E+01
 EBVSHRINKSD(%)  5.4359E+00  6.5574E+01  9.9265E+01  1.4998E+01  3.9706E+01
 EBVSHRINKVR(%)  1.0576E+01  8.8149E+01  9.9995E+01  2.7747E+01  6.3646E+01
 RELATIVEINF(%)  8.6681E+01  5.0897E+00  1.7336E-04  1.4446E+01  1.1284E+00
 EPSSHRINKSD(%)  2.9409E+01
 EPSSHRINKVR(%)  5.0169E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1298.8859099410952     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -563.73508337735700     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1298.886       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.71E-01  2.59E-01  1.22E+00  2.50E-01  8.07E-01  2.78E+00  1.00E-02  9.73E-01  5.65E-01  3.33E+00
 


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
+        1.49E+03
 
 TH 2
+       -2.13E+01  3.67E+03
 
 TH 3
+       -1.15E+02 -1.73E+03  1.14E+04
 
 TH 4
+       -7.63E+01 -4.47E+01 -4.64E+02  5.86E+02
 
 TH 5
+        4.72E+02  1.20E+03 -1.57E+04 -5.14E+02  2.47E+04
 
 TH 6
+       -1.16E+00  8.06E+01 -4.74E+00 -2.46E+01  6.62E+01  2.45E+02
 
 TH 7
+        8.30E+00  3.27E+02 -2.54E+02 -1.52E+01  3.23E+02  5.90E+00  3.18E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.51E+01  7.89E+01  5.39E+01 -1.91E+01  1.08E+02  6.10E+00  5.19E+00  0.00E+00  1.11E+02
 
 TH10
+       -2.33E+01 -2.76E+02 -6.85E+01  1.24E+01  1.45E+02 -5.30E+00 -1.62E+01  0.00E+00 -4.22E+00  1.11E+02
 
 TH11
+       -2.79E+01 -1.95E+02  1.45E+02  1.82E+00 -1.82E+02  1.14E-01 -1.91E+01  0.00E+00  3.78E+00  3.01E+01  4.29E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.413
Stop Time:
Wed Sep 29 13:38:23 CDT 2021
