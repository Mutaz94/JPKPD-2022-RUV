Sat Sep 25 07:37:49 CDT 2021
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
$DATA ../../../../data/spa/B/dat92.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1653.64665841387        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.7441E+01 -7.8358E+01 -6.4574E+01 -3.1336E+01  7.7399E+01  1.9576E+01 -1.9166E+01  9.0652E+00 -1.4773E+01  1.6648E+01
             1.6009E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1664.36671994058        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8331E-01  1.1132E+00  1.1593E+00  9.7545E-01  1.0396E+00  9.6413E-01  1.1426E+00  9.6202E-01  1.0823E+00  8.7229E-01
             1.0187E+00
 PARAMETER:  8.3164E-02  2.0723E-01  2.4782E-01  7.5147E-02  1.3888E-01  6.3467E-02  2.3327E-01  6.1277E-02  1.7909E-01 -3.6639E-02
             1.1852E-01
 GRADIENT:   4.0011E+01  9.0000E+00  3.5638E+00  1.3209E+01  4.7804E+00  6.6139E+00 -2.2650E+00 -2.8421E+00  4.3403E+00 -3.1931E+00
             3.6772E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1664.64894149135        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      205
 NPARAMETR:  9.8290E-01  1.1508E+00  1.1207E+00  9.5222E-01  1.0406E+00  9.5541E-01  1.2040E+00  9.9328E-01  1.0517E+00  8.5623E-01
             1.0106E+00
 PARAMETER:  8.2750E-02  2.4042E-01  2.1394E-01  5.1037E-02  1.3982E-01  5.4389E-02  2.8562E-01  9.3261E-02  1.5044E-01 -5.5217E-02
             1.1057E-01
 GRADIENT:  -1.3672E+00  6.9062E+00  2.3110E+00  7.6017E+00  3.7907E+00 -1.0812E+00  1.5831E+00 -1.1528E+00  3.1846E-03 -3.4326E+00
             4.1149E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.95332261164        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      383
 NPARAMETR:  9.8680E-01  1.3767E+00  8.8100E-01  8.0407E-01  1.0473E+00  9.5970E-01  1.0352E+00  8.2389E-01  1.1684E+00  8.6318E-01
             1.0067E+00
 PARAMETER:  8.6713E-02  4.1972E-01 -2.6699E-02 -1.1807E-01  1.4624E-01  5.8867E-02  1.3462E-01 -9.3713E-02  2.5562E-01 -4.7134E-02
             1.0673E-01
 GRADIENT:   4.1700E+00  1.3436E+01  3.8608E-01  1.1366E+01 -5.5185E+00  1.9346E-03 -2.0197E+00  5.1963E-01 -2.5049E+00  4.1909E-01
            -4.1573E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1665.93047852023        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  9.8194E-01  1.7584E+00  5.8913E-01  5.3286E-01  1.1441E+00  9.6089E-01  8.6978E-01  3.5669E-01  1.5742E+00  8.9779E-01
             1.0112E+00
 PARAMETER:  8.1777E-02  6.6442E-01 -4.2911E-01 -5.2949E-01  2.3458E-01  6.0104E-02 -3.9518E-02 -9.3089E-01  5.5374E-01 -7.8215E-03
             1.1109E-01
 GRADIENT:  -1.0167E+01 -6.3261E+00 -2.5788E+00  1.7119E+00  3.0293E+00  4.4816E-02 -1.0822E+00  2.3723E-01  7.2537E-01  9.9075E-01
             1.3718E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1666.08353428246        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  9.8633E-01  1.8768E+00  5.2758E-01  4.5384E-01  1.1877E+00  9.6073E-01  8.3978E-01  2.4507E-01  1.7477E+00  9.1560E-01
             1.0099E+00
 PARAMETER:  8.6239E-02  7.2957E-01 -5.3946E-01 -6.9000E-01  2.7200E-01  5.9940E-02 -7.4610E-02 -1.3062E+00  6.5828E-01  1.1826E-02
             1.0986E-01
 GRADIENT:   3.6038E-01 -1.9243E+00 -3.6742E-01 -6.9220E-01  1.6909E-01  3.5096E-03  2.9986E-01  1.0705E-01  2.1669E-02  1.3008E-01
             7.2532E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1666.11519525985        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      914
 NPARAMETR:  9.8585E-01  1.8621E+00  5.4096E-01  4.6504E-01  1.1817E+00  9.6087E-01  8.4657E-01  1.4234E-01  1.7339E+00  9.1566E-01
             1.0100E+00
 PARAMETER:  8.5745E-02  7.2173E-01 -5.1440E-01 -6.6564E-01  2.6697E-01  6.0080E-02 -6.6561E-02 -1.8495E+00  6.5036E-01  1.1895E-02
             1.0994E-01
 GRADIENT:  -7.3296E-01  1.8943E+00  9.6024E-01 -4.4291E-01 -3.0665E+00  6.7376E-02  8.7434E-01  3.2371E-02  6.1351E-01  3.7981E-01
             8.1141E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1666.14069468552        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1089
 NPARAMETR:  9.8628E-01  1.8578E+00  5.3657E-01  4.6736E-01  1.1792E+00  9.6076E-01  8.4543E-01  1.9960E-02  1.7175E+00  9.1136E-01
             1.0097E+00
 PARAMETER:  8.6190E-02  7.1941E-01 -5.2255E-01 -6.6066E-01  2.6483E-01  5.9966E-02 -6.7914E-02 -3.8140E+00  6.4087E-01  7.1808E-03
             1.0969E-01
 GRADIENT:   2.4592E-01 -1.0231E-01  2.2346E-02 -4.2602E-02 -9.9323E-02  1.2090E-02  5.9431E-02  6.9178E-04  3.7029E-02  1.1275E-02
            -3.6905E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1666.14094050803        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1271             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8619E-01  1.8590E+00  5.3586E-01  4.6658E-01  1.1797E+00  9.6073E-01  8.4477E-01  1.0000E-02  1.7191E+00  9.1159E-01
             1.0098E+00
 PARAMETER:  8.6096E-02  7.2005E-01 -5.2389E-01 -6.6232E-01  2.6530E-01  5.9940E-02 -6.8695E-02 -4.7107E+00  6.4181E-01  7.4343E-03
             1.0980E-01
 GRADIENT:   4.1096E+01  8.0031E+01  2.3823E-01  1.0215E+01  1.6552E+00  4.1534E+00  8.5822E-01  0.0000E+00  2.6951E+00  5.3538E-02
             8.8735E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1666.14094194572        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1364
 NPARAMETR:  9.8619E-01  1.8591E+00  5.3589E-01  4.6658E-01  1.1797E+00  9.6073E-01  8.4485E-01  1.0000E-02  1.7191E+00  9.1160E-01
             1.0098E+00
 PARAMETER:  8.6092E-02  7.2009E-01 -5.2383E-01 -6.6232E-01  2.6527E-01  5.9937E-02 -6.8595E-02 -4.7107E+00  6.4179E-01  7.4420E-03
             1.0980E-01
 GRADIENT:   6.9736E-03 -2.4669E-02  7.2059E-03 -6.3214E-03 -6.9385E-03  1.3246E-03  1.3779E-03  0.0000E+00  1.3168E-03  9.1676E-04
            -7.0423E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1364
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1483E-04 -3.0738E-02 -2.3730E-04  3.1117E-02 -3.9780E-02
 SE:             2.9850E-02  2.4597E-02  8.8039E-05  2.2144E-02  2.1714E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9693E-01  2.1143E-01  7.0296E-03  1.5997E-01  6.6948E-02

 ETASHRINKSD(%)  1.0000E-10  1.7595E+01  9.9705E+01  2.5814E+01  2.7257E+01
 ETASHRINKVR(%)  1.0000E-10  3.2095E+01  9.9999E+01  4.4965E+01  4.7085E+01
 EBVSHRINKSD(%)  4.6924E-01  1.6603E+01  9.9749E+01  2.8374E+01  2.5864E+01
 EBVSHRINKVR(%)  9.3627E-01  3.0450E+01  9.9999E+01  4.8697E+01  4.5038E+01
 RELATIVEINF(%)  9.9018E+01  6.9870E+00  9.9895E-05  4.9243E+00  1.4953E+01
 EPSSHRINKSD(%)  4.4174E+01
 EPSSHRINKVR(%)  6.8834E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1666.1409419457161     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.99011538197794     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1666.141       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.86E+00  5.36E-01  4.67E-01  1.18E+00  9.61E-01  8.45E-01  1.00E-02  1.72E+00  9.12E-01  1.01E+00
 


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
+        1.23E+03
 
 TH 2
+       -5.33E+00  3.30E+02
 
 TH 3
+        7.55E+00  1.17E+02  2.83E+02
 
 TH 4
+       -1.46E+01  2.66E+02 -2.67E+02  9.29E+02
 
 TH 5
+       -5.37E+00 -1.49E+02 -2.46E+02  2.70E+02  5.31E+02
 
 TH 6
+       -7.26E-01 -5.48E-01  1.90E+00 -2.03E+00  9.65E-02  2.16E+02
 
 TH 7
+        1.38E+00  8.30E+00 -1.07E+01 -1.52E+01 -1.42E+01 -2.90E+00  1.47E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.88E-01 -1.48E+01 -2.93E+01  5.90E+01 -1.66E-01 -3.98E-01  1.41E+01  0.00E+00  2.59E+01
 
 TH10
+        2.56E+00 -1.42E+01 -2.54E+01 -3.99E+00 -6.70E+01  2.07E-01  9.07E+00  0.00E+00  5.22E+00  9.67E+01
 
 TH11
+       -5.81E+00 -1.55E+01 -1.84E+01  1.04E+00 -8.54E+00  3.43E+00  9.09E+00  0.00E+00  3.94E+00  1.51E+01  2.13E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.473
Stop Time:
Sat Sep 25 07:38:14 CDT 2021
