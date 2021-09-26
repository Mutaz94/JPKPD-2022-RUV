Sat Sep 25 10:15:59 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat2.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1672.23335482836        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.6488E+01 -6.0480E+01 -2.6831E+01 -6.0992E+01 -9.2760E+00  7.5624E+00  1.0682E+01  1.5491E+01  3.1081E+01  1.7574E+01
            -3.3306E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1676.78088044913        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8620E-01  1.1564E+00  1.3019E+00  9.6735E-01  1.1991E+00  9.7321E-01  8.8836E-01  8.3141E-01  7.2124E-01  9.0999E-01
             1.1833E+00
 PARAMETER:  8.6103E-02  2.4528E-01  3.6382E-01  6.6800E-02  2.8155E-01  7.2844E-02 -1.8381E-02 -8.4627E-02 -2.2678E-01  5.6815E-03
             2.6829E-01
 GRADIENT:   4.2421E+01  2.0918E+01  1.5526E+01  5.0782E+00  1.7297E+01 -1.6731E+00 -1.2966E+01 -4.3206E+00 -2.5600E+01 -2.7199E+01
             9.1382E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1681.93011611835        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8712E-01  9.9457E-01  1.1341E+00  1.0620E+00  1.0910E+00  1.0155E+00  8.7208E-01  2.1288E-01  9.0843E-01  1.0442E+00
             1.1601E+00
 PARAMETER:  8.7032E-02  9.4559E-02  2.2583E-01  1.6020E-01  1.8712E-01  1.1542E-01 -3.6876E-02 -1.4470E+00  3.9645E-03  1.4330E-01
             2.4852E-01
 GRADIENT:   4.5431E+01 -1.1713E+01 -1.3505E+01  1.6786E+01  3.3308E+01  1.4610E+01  5.4853E-01 -1.3568E-01  1.1718E+01 -1.5886E+00
             2.1720E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1683.84948664193        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  9.6700E-01  1.0600E+00  9.6933E-01  1.0084E+00  1.0125E+00  9.8043E-01  9.7568E-01  2.8026E-01  8.3697E-01  9.5139E-01
             1.0935E+00
 PARAMETER:  6.6441E-02  1.5825E-01  6.8850E-02  1.0837E-01  1.1244E-01  8.0236E-02  7.5374E-02 -1.1720E+00 -7.7966E-02  5.0172E-02
             1.8936E-01
 GRADIENT:  -3.2205E+01 -3.2303E+00 -3.4759E-01 -5.1154E+00  2.4709E-01 -2.8909E+00  2.1544E-01  4.1068E-02 -5.5012E-01 -2.8717E-01
             3.0668E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1684.10810871257        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      430
 NPARAMETR:  9.8092E-01  1.1112E+00  9.6029E-01  9.8059E-01  1.0318E+00  9.8631E-01  9.3056E-01  2.8369E-01  8.6545E-01  9.6376E-01
             1.0936E+00
 PARAMETER:  8.0732E-02  2.0545E-01  5.9475E-02  8.0396E-02  1.3132E-01  8.6217E-02  2.8035E-02 -1.1599E+00 -4.4504E-02  6.3087E-02
             1.8944E-01
 GRADIENT:  -1.5429E-01  1.7615E+00  4.8204E-01  1.8445E+00 -1.1622E+00 -2.1820E-01  1.0425E-01  2.2136E-02  1.3274E-01 -4.7062E-02
            -6.7707E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1684.11535047854        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  9.8135E-01  1.1744E+00  9.3426E-01  9.3860E-01  1.0524E+00  9.8735E-01  8.8886E-01  2.6937E-01  8.9620E-01  9.7150E-01
             1.0938E+00
 PARAMETER:  8.1178E-02  2.6079E-01  3.2002E-02  3.6634E-02  1.5107E-01  8.7268E-02 -1.7814E-02 -1.2117E+00 -9.5926E-03  7.1088E-02
             1.8969E-01
 GRADIENT:   7.6023E-02 -6.6692E-01 -1.6866E-01 -6.0619E-01  1.8456E-01  5.3643E-02  7.4141E-02  3.7316E-02  1.7429E-01  9.4868E-02
             1.1332E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1684.12597271243        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      782
 NPARAMETR:  9.8177E-01  1.1948E+00  9.0469E-01  9.2485E-01  1.0473E+00  9.8821E-01  8.8239E-01  1.1818E-01  9.0109E-01  9.6589E-01
             1.0943E+00
 PARAMETER:  8.1602E-02  2.7796E-01 -1.6466E-04  2.1874E-02  1.4623E-01  8.8139E-02 -2.5116E-02 -2.0355E+00 -4.1544E-03  6.5290E-02
             1.9015E-01
 GRADIENT:   4.2712E-01  9.8023E-02 -1.6545E-01  6.8360E-02 -2.1513E-01  2.8767E-01 -5.0903E-02  6.6344E-03 -1.5992E-01  1.9246E-01
             2.4982E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1684.12839552759        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  9.8162E-01  1.2016E+00  9.0050E-01  9.2027E-01  1.0490E+00  9.8754E-01  8.7901E-01  3.8107E-02  9.0528E-01  9.6573E-01
             1.0940E+00
 PARAMETER:  8.1446E-02  2.8361E-01 -4.8072E-03  1.6912E-02  1.4784E-01  8.7458E-02 -2.8964E-02 -3.1674E+00  4.8952E-04  6.5133E-02
             1.8986E-01
 GRADIENT:  -1.4089E-03 -5.3632E-02  4.2901E-02 -9.7848E-02  4.8870E-02  2.0487E-03 -7.1190E-03  3.7683E-04 -2.4210E-02 -4.0997E-02
            -2.8048E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1684.12861566638        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1134
 NPARAMETR:  9.8162E-01  1.2000E+00  9.0017E-01  9.2128E-01  1.0479E+00  9.8753E-01  8.8014E-01  1.0000E-02  9.0446E-01  9.6528E-01
             1.0940E+00
 PARAMETER:  8.1447E-02  2.8229E-01 -5.1710E-03  1.8010E-02  1.4680E-01  8.7454E-02 -2.7670E-02 -4.7524E+00 -4.1470E-04  6.4664E-02
             1.8986E-01
 GRADIENT:   1.5005E-03 -7.1502E-03 -1.0760E-03 -7.7450E-03 -3.1134E-04  3.7161E-04  1.7434E-03  0.0000E+00  2.4307E-03  1.9856E-03
             3.0134E-04

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1684.12861566638        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1156
 NPARAMETR:  9.8162E-01  1.2000E+00  9.0017E-01  9.2128E-01  1.0479E+00  9.8753E-01  8.8014E-01  1.0000E-02  9.0446E-01  9.6528E-01
             1.0940E+00
 PARAMETER:  8.1447E-02  2.8229E-01 -5.1710E-03  1.8010E-02  1.4680E-01  8.7454E-02 -2.7670E-02 -4.7524E+00 -4.1470E-04  6.4664E-02
             1.8986E-01
 GRADIENT:   1.5005E-03 -7.1502E-03 -1.0760E-03 -7.7450E-03 -3.1134E-04  3.7161E-04  1.7434E-03  0.0000E+00  2.4307E-03  1.9856E-03
             3.0134E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1156
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.5810E-04 -1.2988E-02 -2.9423E-04  2.5017E-03 -2.6186E-02
 SE:             2.9781E-02  2.0006E-02  1.4307E-04  2.3446E-02  2.3137E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9308E-01  5.1620E-01  3.9733E-02  9.1503E-01  2.5773E-01

 ETASHRINKSD(%)  2.2984E-01  3.2979E+01  9.9521E+01  2.1453E+01  2.2488E+01
 ETASHRINKVR(%)  4.5915E-01  5.5081E+01  9.9998E+01  3.8304E+01  3.9920E+01
 EBVSHRINKSD(%)  5.1768E-01  3.2350E+01  9.9553E+01  2.2082E+01  2.0935E+01
 EBVSHRINKVR(%)  1.0327E+00  5.4235E+01  9.9998E+01  3.9288E+01  3.7487E+01
 RELATIVEINF(%)  9.8467E+01  1.1566E+00  1.7664E-04  1.7759E+00  6.5482E+00
 EPSSHRINKSD(%)  4.1412E+01
 EPSSHRINKVR(%)  6.5675E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1684.1286156663841     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -948.97778910264594     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.14
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1684.129       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.20E+00  9.00E-01  9.21E-01  1.05E+00  9.88E-01  8.80E-01  1.00E-02  9.04E-01  9.65E-01  1.09E+00
 


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
+        1.17E+03
 
 TH 2
+       -1.00E+01  4.73E+02
 
 TH 3
+        1.50E+01  1.49E+02  2.97E+02
 
 TH 4
+       -1.66E+01  4.80E+02 -1.64E+02  9.76E+02
 
 TH 5
+       -3.61E+00 -2.59E+02 -3.73E+02  1.79E+02  6.77E+02
 
 TH 6
+        1.13E+00 -2.17E+00  6.85E+00 -5.30E+00 -1.30E+00  1.99E+02
 
 TH 7
+        2.75E+00  1.80E+01  7.61E+00 -1.18E+01 -1.37E+01  1.58E+00  5.49E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.23E+00 -1.89E+01 -1.91E+01  2.52E+01  6.62E+00  3.21E+00  3.41E+01  0.00E+00  9.55E+01
 
 TH10
+       -2.99E+00 -5.72E+00 -2.88E+01 -1.34E+01 -5.53E+01  1.32E+00  1.28E+01  0.00E+00  3.40E+00  8.80E+01
 
 TH11
+       -1.04E+01 -2.12E+01 -3.43E+01  4.37E-01  6.34E+00  2.13E+00  7.91E+00  0.00E+00  1.39E+01  2.23E+01  1.91E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.005
Stop Time:
Sat Sep 25 10:16:19 CDT 2021
