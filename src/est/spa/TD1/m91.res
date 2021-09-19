Sat Sep 18 14:21:05 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat91.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1606.92297392677        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8483E+02 -8.3598E+01 -5.8707E+01 -6.5846E+01  4.6675E+01 -1.8782E+01 -4.7328E+00  1.8603E+01 -8.9156E+00  2.7384E+01
            -4.7671E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1617.23505745251        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9276E-01  1.0563E+00  1.2771E+00  9.8386E-01  1.1172E+00  1.0873E+00  1.0078E+00  8.4642E-01  1.0772E+00  7.9801E-01
             1.2088E+00
 PARAMETER:  9.2732E-02  1.5473E-01  3.4457E-01  8.3730E-02  2.1086E-01  1.8367E-01  1.0780E-01 -6.6741E-02  1.7441E-01 -1.2564E-01
             2.8959E-01
 GRADIENT:   1.3618E+02 -5.1833E+01 -8.1062E+00 -5.3523E+01  6.3929E+01  1.8492E+01  2.7303E+00  1.9286E+00  1.8875E+00 -1.7780E+01
             2.4299E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1621.97662148786        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.8337E-01  9.0153E-01  1.1654E+00  1.1159E+00  1.0143E+00  1.0671E+00  7.9583E-01  2.7268E-01  9.9968E-01  9.5443E-01
             1.1794E+00
 PARAMETER:  8.3229E-02 -3.6670E-03  2.5306E-01  2.0970E-01  1.1423E-01  1.6497E-01 -1.2838E-01 -1.1995E+00  9.9678E-02  5.3364E-02
             2.6504E-01
 GRADIENT:   1.2277E+02 -2.1483E+01 -2.6917E+01 -1.2698E+00  5.1884E+01  1.3664E+01 -4.1769E+00 -1.4637E-01 -1.1632E+01  1.7243E+00
             2.2453E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1627.79096076661        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      274
 NPARAMETR:  9.3346E-01  8.3381E-01  1.0908E+00  1.1636E+00  9.3696E-01  1.0100E+00  1.2558E+00  2.6312E-01  9.5734E-01  8.7923E-01
             1.1023E+00
 PARAMETER:  3.1138E-02 -8.1745E-02  1.8693E-01  2.5149E-01  3.4887E-02  1.0996E-01  3.2780E-01 -1.2351E+00  5.6402E-02 -2.8711E-02
             1.9738E-01
 GRADIENT:  -3.6885E+00 -5.5043E+00 -9.6904E+00  5.9708E+00  2.1800E+01 -7.5181E+00  5.5108E-01  2.1442E-01  6.1656E+00  6.9544E+00
            -2.2476E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1629.64637863630        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.3702E-01  6.9924E-01  9.1511E-01  1.2273E+00  7.8599E-01  1.0270E+00  1.5771E+00  1.5681E-01  8.3907E-01  6.9329E-01
             1.1140E+00
 PARAMETER:  3.4944E-02 -2.5777E-01  1.1293E-02  3.0483E-01 -1.4081E-01  1.2669E-01  5.5559E-01 -1.7527E+00 -7.5459E-02 -2.6631E-01
             2.0793E-01
 GRADIENT:   3.6791E+00  4.3721E+00  1.0922E+00  4.2273E-01 -3.2321E+00 -9.6713E-01  4.4654E-01  3.1526E-01 -3.4707E+00  4.7548E-01
             1.9031E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1630.08109039184        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      626
 NPARAMETR:  9.3446E-01  5.6880E-01  9.0834E-01  1.2992E+00  7.4165E-01  1.0270E+00  1.8478E+00  2.1859E-02  8.1389E-01  6.8773E-01
             1.1073E+00
 PARAMETER:  3.2209E-02 -4.6423E-01  3.8671E-03  3.6177E-01 -1.9887E-01  1.2662E-01  7.1400E-01 -3.7231E+00 -1.0593E-01 -2.7435E-01
             2.0194E-01
 GRADIENT:   7.8300E-01  1.2480E+00  1.6879E+00 -1.3736E-01 -3.0712E+00 -4.7339E-01  3.5723E-01  5.4542E-03  4.0915E-01  2.1023E-01
             8.4554E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1630.09636170215        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      802            RESET HESSIAN, TYPE II
 NPARAMETR:  9.3375E-01  5.3849E-01  9.1068E-01  1.3160E+00  7.3489E-01  1.0278E+00  1.9212E+00  1.0000E-02  8.0401E-01  6.8890E-01
             1.1070E+00
 PARAMETER:  3.1459E-02 -5.1899E-01  6.4367E-03  3.7461E-01 -2.0803E-01  1.2743E-01  7.5293E-01 -4.5734E+00 -1.1814E-01 -2.7266E-01
             2.0163E-01
 GRADIENT:   2.9711E+01  4.0222E+00  3.1456E-01  3.1751E+01  1.6542E+00  3.8576E+00  1.9874E+00  0.0000E+00  6.9469E-01  1.0108E-01
             1.5913E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1630.09636275158        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.3372E-01  5.3855E-01  9.1066E-01  1.3159E+00  7.3490E-01  1.0278E+00  1.9211E+00  1.0000E-02  8.0401E-01  6.8887E-01
             1.1070E+00
 PARAMETER:  3.1426E-02 -5.1888E-01  6.4109E-03  3.7455E-01 -2.0803E-01  1.2738E-01  7.5290E-01 -5.0744E+00 -1.1814E-01 -2.7270E-01
             2.0163E-01
 GRADIENT:   2.9640E+01  4.0131E+00  3.2769E-01  3.1662E+01  1.6419E+00  3.8377E+00  1.9937E+00  0.0000E+00  6.9052E-01  1.0179E-01
             1.6142E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1630.09639498561        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1073
 NPARAMETR:  9.3376E-01  5.3948E-01  9.1034E-01  1.3154E+00  7.3500E-01  1.0278E+00  1.9182E+00  1.0000E-02  8.0431E-01  6.8868E-01
             1.1069E+00
 PARAMETER:  3.1466E-02 -5.1715E-01  6.0656E-03  3.7417E-01 -2.0788E-01  1.2747E-01  7.5136E-01 -6.3025E+00 -1.1777E-01 -2.7298E-01
             2.0161E-01
 GRADIENT:  -8.9063E-03 -8.2256E-03 -8.0211E-03  3.2225E-02 -1.7951E-02  9.6175E-03 -1.6317E-02  0.0000E+00  2.1763E-03 -1.0980E-03
            -1.2104E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1630.09640205827        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1174
 NPARAMETR:  9.3380E-01  5.3979E-01  9.1026E-01  1.3152E+00  7.3509E-01  1.0279E+00  1.9166E+00  1.0000E-02  8.0429E-01  6.8848E-01
             1.1069E+00
 PARAMETER:  3.1474E-02 -5.1645E-01  5.9448E-03  3.7397E-01 -2.0781E-01  1.2743E-01  7.5111E-01 -6.9617E+00 -1.1768E-01 -2.7312E-01
             2.0165E-01
 GRADIENT:  -1.2561E-02  5.5742E-03 -7.8842E-03  1.7249E-03 -2.2430E-02 -8.2256E-03  1.1358E-02  0.0000E+00  3.8563E-03  3.1264E-03
             8.0286E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1174
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3898E-04  2.2485E-02 -3.9388E-04 -1.9920E-02 -1.8165E-03
 SE:             2.9764E-02  1.9747E-02  2.1170E-04  2.5120E-02  2.1602E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8823E-01  2.5486E-01  6.2807E-02  4.2778E-01  9.3299E-01

 ETASHRINKSD(%)  2.8843E-01  3.3845E+01  9.9291E+01  1.5843E+01  2.7631E+01
 ETASHRINKVR(%)  5.7602E-01  5.6235E+01  9.9995E+01  2.9177E+01  4.7627E+01
 EBVSHRINKSD(%)  5.3959E-01  3.6302E+01  9.9259E+01  1.4633E+01  2.5331E+01
 EBVSHRINKVR(%)  1.0763E+00  5.9425E+01  9.9995E+01  2.7125E+01  4.4245E+01
 RELATIVEINF(%)  9.8433E+01  4.9502E+00  4.7123E-04  1.2278E+01  3.9338E+00
 EPSSHRINKSD(%)  4.1432E+01
 EPSSHRINKVR(%)  6.5698E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1630.0964020582724     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -894.94557549453418     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1630.096       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.34E-01  5.40E-01  9.10E-01  1.32E+00  7.35E-01  1.03E+00  1.92E+00  1.00E-02  8.04E-01  6.89E-01  1.11E+00
 


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
+        1.19E+03
 
 TH 2
+       -1.50E+01  4.56E+02
 
 TH 3
+        1.58E+01  2.65E+02  7.78E+02
 
 TH 4
+       -7.24E+00  3.10E+02 -2.14E+02  7.38E+02
 
 TH 5
+       -5.64E+00 -5.69E+02 -1.27E+03  2.67E+02  2.40E+03
 
 TH 6
+        9.33E-01 -3.39E+00  5.97E+00 -3.54E+00 -6.78E+00  1.82E+02
 
 TH 7
+        1.05E+00  3.74E+01  3.02E+00 -7.15E+00 -5.65E+00  4.72E-02  1.59E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.38E+00 -2.39E+01 -2.31E+01  3.76E+00  1.75E+01 -1.49E+00  1.03E+01  0.00E+00  1.69E+02
 
 TH10
+       -2.98E+00  5.48E+00 -6.49E+01 -3.50E+01 -3.44E+01 -8.89E-01  6.20E+00  0.00E+00  7.10E+00  1.33E+02
 
 TH11
+       -9.40E+00 -9.10E+00 -4.18E+01 -8.29E+00  2.04E+01  2.70E+00  2.25E+00  0.00E+00  1.03E+01  3.58E+01  1.86E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.981
Stop Time:
Sat Sep 18 14:21:26 CDT 2021
