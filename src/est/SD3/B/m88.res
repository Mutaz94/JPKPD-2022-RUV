Sat Oct 23 21:25:58 CDT 2021
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
$DATA ../../../../data/SD3/B/dat88.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2134.87495820080        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1759E+02 -4.6496E+01 -4.3476E+01  1.1293E+01  1.0586E+02  4.9796E+01 -3.9888E+00 -1.0075E+01 -5.2394E+00  1.3723E+01
             4.2749E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2143.25655900518        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  1.0154E+00  1.0317E+00  1.0201E+00  1.0533E+00  9.5785E-01  9.7317E-01  1.0166E+00  1.0892E+00  1.0571E+00  8.8008E-01
             9.1225E-01
 PARAMETER:  1.1528E-01  1.3116E-01  1.1990E-01  1.5189E-01  5.6940E-02  7.2803E-02  1.1650E-01  1.8541E-01  1.5552E-01 -2.7737E-02
             8.1561E-03
 GRADIENT:   2.7266E+01  8.6582E-01 -2.1576E+01  3.3271E+01  4.2418E+01 -4.4913E+00 -6.8208E-01 -9.0908E+00  4.1377E+00  2.8693E+00
            -3.3911E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2146.00430599952        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  1.0154E+00  9.5943E-01  1.2300E+00  1.1029E+00  9.9342E-01  9.8235E-01  1.1127E+00  1.5150E+00  1.0222E+00  8.2903E-01
             9.3446E-01
 PARAMETER:  1.1531E-01  5.8582E-02  3.0699E-01  1.9798E-01  9.3398E-02  8.2188E-02  2.0675E-01  5.1543E-01  1.2198E-01 -8.7497E-02
             3.2216E-02
 GRADIENT:   2.8756E+01  1.5114E+00 -1.8218E+01  1.8915E+01  2.6905E+01 -3.4400E-01  4.5165E+00  5.1464E+00  6.0490E+00 -2.1247E+00
            -1.6011E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2148.40550213879        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.9886E-01  7.8500E-01  1.7633E+00  1.2148E+00  1.0558E+00  9.8183E-01  7.3334E-01  1.8679E+00  1.0084E+00  9.3246E-01
             9.5655E-01
 PARAMETER:  9.8862E-02 -1.4207E-01  6.6721E-01  2.9461E-01  1.5432E-01  8.1666E-02 -2.1014E-01  7.2479E-01  1.0836E-01  3.0075E-02
             5.5581E-02
 GRADIENT:  -5.8102E+00  3.6405E+00  3.4742E+00 -4.9661E-01 -8.1149E+00  6.1044E-01  1.9171E+00  8.2710E-01  2.9562E+00  2.4385E+00
             2.9573E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2149.42223488047        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.9945E-01  5.6512E-01  1.8293E+00  1.3577E+00  1.0139E+00  9.7724E-01  1.6615E-01  1.8193E+00  9.3508E-01  9.0472E-01
             9.5077E-01
 PARAMETER:  9.9447E-02 -4.7072E-01  7.0392E-01  4.0583E-01  1.1385E-01  7.6977E-02 -1.6948E+00  6.9845E-01  3.2875E-02 -1.2622E-04
             4.9521E-02
 GRADIENT:   8.3590E-01  2.8792E+00 -1.5073E-01  7.7561E+00 -3.7416E-01 -2.2293E-01  7.2020E-02 -3.4097E-01 -8.2233E-02 -3.6374E-01
            -7.3245E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2149.46931577043        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.9873E-01  5.5667E-01  1.8264E+00  1.3584E+00  1.0113E+00  9.7765E-01  5.5246E-02  1.8156E+00  9.3293E-01  9.0609E-01
             9.5132E-01
 PARAMETER:  9.8727E-02 -4.8579E-01  7.0237E-01  4.0630E-01  1.1127E-01  7.7399E-02 -2.7960E+00  6.9641E-01  3.0580E-02  1.3868E-03
             5.0093E-02
 GRADIENT:  -5.3381E-01  4.1880E-01 -2.7281E-01 -4.4841E-01 -2.2685E-01  2.4321E-03  1.0312E-02 -1.1565E-01  1.5987E-01  8.4954E-03
             3.0581E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2149.47659319358        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1072             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0000E+00  5.5619E-01  1.8323E+00  1.3575E+00  1.0118E+00  9.7792E-01  1.0000E-02  1.8217E+00  9.3171E-01  9.0588E-01
             9.5130E-01
 PARAMETER:  1.0001E-01 -4.8664E-01  7.0555E-01  4.0564E-01  1.1177E-01  7.7669E-02 -5.8305E+00  6.9978E-01  2.9266E-02  1.1566E-03
             5.0079E-02
 GRADIENT:   4.7688E+02  7.1479E+01  1.1986E+01  5.7660E+02  7.0929E+00  4.6784E+01  0.0000E+00  6.8746E+00  1.1447E+01  5.9651E-01
             9.4124E-01

0ITERATION NO.:   34    OBJECTIVE VALUE:  -2149.47748827216        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1206
 NPARAMETR:  9.9974E-01  5.5628E-01  1.8318E+00  1.3567E+00  1.0124E+00  9.7790E-01  1.0000E-02  1.8222E+00  9.3265E-01  9.0591E-01
             9.5137E-01
 PARAMETER:  9.9735E-02 -4.8648E-01  7.0529E-01  4.0504E-01  1.1232E-01  7.7651E-02 -5.8305E+00  7.0007E-01  3.0276E-02  1.1902E-03
             5.0152E-02
 GRADIENT:   1.8117E+00 -6.6873E-01 -2.8139E-01 -4.3339E+00  1.1163E-01  1.1109E-01  0.0000E+00  7.5832E-02  9.7113E-02 -6.6788E-02
             1.7550E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1206
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.7160E-04 -2.8230E-04 -4.2696E-02 -4.4805E-03 -4.3542E-02
 SE:             2.9898E-02  1.0424E-04  1.9426E-02  2.9613E-02  2.0289E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9008E-01  6.7662E-03  2.7956E-02  8.7974E-01  3.1867E-02

 ETASHRINKSD(%)  1.0000E-10  9.9651E+01  3.4921E+01  7.9319E-01  3.2028E+01
 ETASHRINKVR(%)  1.0000E-10  9.9999E+01  5.7647E+01  1.5801E+00  5.3798E+01
 EBVSHRINKSD(%)  3.0912E-01  9.9681E+01  3.7502E+01  1.4048E+00  2.9349E+01
 EBVSHRINKVR(%)  6.1728E-01  9.9999E+01  6.0940E+01  2.7899E+00  5.0084E+01
 RELATIVEINF(%)  9.8641E+01  8.1994E-05  1.3065E+01  9.0488E+00  1.4118E+01
 EPSSHRINKSD(%)  3.4733E+01
 EPSSHRINKVR(%)  5.7403E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2149.4774882721604     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1230.5389550674877     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2149.477       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  5.56E-01  1.83E+00  1.36E+00  1.01E+00  9.78E-01  1.00E-02  1.82E+00  9.33E-01  9.06E-01  9.51E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.01
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       98.325
Stop Time:
Sat Oct 23 21:26:14 CDT 2021
