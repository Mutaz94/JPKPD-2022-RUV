Sat Oct 23 21:59:45 CDT 2021
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
$DATA ../../../../data/SD3/A1/dat93.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m93.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1889.81106792371        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4266E+02 -7.3289E+00  4.2251E+01  2.0303E+01  6.5815E+01  4.2652E+01  1.7564E+01 -8.2707E+01  4.5653E+01 -1.8139E+01
            -3.5149E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1953.28494803443        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8509E-01  1.0583E+00  8.1185E-01  9.9195E-01  8.7490E-01  9.7543E-01  8.5677E-01  1.3394E+00  7.2968E-01  9.1499E-01
             1.6926E+00
 PARAMETER:  8.4975E-02  1.5668E-01 -1.0844E-01  9.1921E-02 -3.3651E-02  7.5119E-02 -5.4581E-02  3.9224E-01 -2.1514E-01  1.1161E-02
             6.2629E-01
 GRADIENT:   1.3264E+02  4.7641E+00  4.0491E+00  4.1118E+00 -1.6371E+01  1.5140E+01 -2.3011E+00  4.9533E+00 -1.0041E+01  2.2337E+01
             1.3091E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1959.59497819616        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  9.8158E-01  1.1435E+00  5.5850E-01  9.3927E-01  7.3345E-01  1.0038E+00  8.8774E-01  1.8108E+00  8.1213E-01  5.8814E-01
             1.5789E+00
 PARAMETER:  8.1411E-02  2.3412E-01 -4.8251E-01  3.7346E-02 -2.1000E-01  1.0375E-01 -1.9074E-02  6.9375E-01 -1.0809E-01 -4.3079E-01
             5.5675E-01
 GRADIENT:  -1.5727E+01  3.7732E+01  2.7212E+00  2.8253E+01 -8.7489E+01  1.0386E+01  4.4256E+00  3.2263E+01  1.1621E+00  1.0446E+01
             1.0394E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1978.13510202421        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      358
 NPARAMETR:  9.8462E-01  1.3899E+00  6.6721E-01  8.0116E-01  9.6503E-01  9.7959E-01  7.7682E-01  1.7859E+00  9.3257E-01  7.5307E-01
             1.3490E+00
 PARAMETER:  8.4499E-02  4.2924E-01 -3.0465E-01 -1.2169E-01  6.4403E-02  7.9378E-02 -1.5254E-01  6.7993E-01  3.0192E-02 -1.8360E-01
             3.9934E-01
 GRADIENT:  -3.5080E+00  1.4799E+01  2.8680E-01  2.6168E+01 -1.1881E+00  6.6824E-01  2.8824E+00 -2.5804E+00  5.8488E-01  1.1008E-01
            -9.1352E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1980.89931547142        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  9.9262E-01  1.6780E+00  6.3021E-01  6.2503E-01  1.1378E+00  9.8576E-01  6.0805E-01  2.0584E+00  1.1729E+00  9.2529E-01
             1.3770E+00
 PARAMETER:  9.2592E-02  6.1763E-01 -3.6170E-01 -3.6996E-01  2.2906E-01  8.5653E-02 -3.9750E-01  8.2194E-01  2.5946E-01  2.2351E-02
             4.1990E-01
 GRADIENT:   1.3494E+01  2.8033E+01  1.1616E-01  2.6204E+01  4.1540E+00  3.4430E+00 -6.2514E+00 -1.0790E+00 -3.1825E+00 -1.4321E+00
             8.8540E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1981.07532991137        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  9.8758E-01  1.6778E+00  6.2572E-01  6.2220E-01  1.1409E+00  9.7809E-01  6.0940E-01  2.0685E+00  1.1793E+00  9.2675E-01
             1.3799E+00
 PARAMETER:  8.7504E-02  6.1750E-01 -3.6886E-01 -3.7449E-01  2.3179E-01  7.7847E-02 -3.9527E-01  8.2683E-01  2.6492E-01  2.3928E-02
             4.2200E-01
 GRADIENT:   2.0608E+00  1.7898E+01 -1.5605E+00  2.3660E+01  8.1993E+00  5.3226E-01 -5.5373E+00 -3.9416E-01 -2.3189E+00 -1.3939E+00
             2.9377E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1981.35643901007        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  9.8277E-01  1.6674E+00  6.2734E-01  6.1535E-01  1.1376E+00  9.7466E-01  6.1913E-01  2.0672E+00  1.1836E+00  9.3077E-01
             1.3774E+00
 PARAMETER:  8.2616E-02  6.1125E-01 -3.6627E-01 -3.8556E-01  2.2894E-01  7.4332E-02 -3.7945E-01  8.2619E-01  2.6855E-01  2.8253E-02
             4.2022E-01
 GRADIENT:  -8.6649E+00 -1.1705E+01 -3.2276E-01  5.7067E+00  5.6108E+00 -7.8385E-01 -3.2490E+00 -1.9639E-01 -7.8508E-01 -5.3182E-01
             2.9844E+00

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1981.53069158400        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  9.8694E-01  1.6738E+00  6.2192E-01  6.1065E-01  1.1362E+00  9.7690E-01  6.3900E-01  2.0798E+00  1.1812E+00  9.3194E-01
             1.3730E+00
 PARAMETER:  8.7004E-02  6.1651E-01 -3.7409E-01 -3.9414E-01  2.2714E-01  7.6739E-02 -3.4865E-01  8.3039E-01  2.6919E-01  2.9736E-02
             4.1602E-01
 GRADIENT:   7.5680E-01  3.7675E+03  3.1218E+03 -5.9186E+03 -1.0264E+04  9.6798E-02 -6.6823E+03 -2.8562E+03  4.0288E-01  1.1657E+04
            -5.6031E+03
 NUMSIGDIG:         2.5         2.3         2.3         2.3         2.3         2.6         2.3         2.3         1.7         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      944
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.5586E-04 -3.0442E-02 -3.2689E-02  2.1642E-02 -4.4708E-02
 SE:             2.9810E-02  2.1082E-02  1.5744E-02  2.2702E-02  1.9363E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8245E-01  1.4873E-01  3.7869E-02  3.4043E-01  2.0947E-02

 ETASHRINKSD(%)  1.3380E-01  2.9374E+01  4.7256E+01  2.3947E+01  3.5130E+01
 ETASHRINKVR(%)  2.6742E-01  5.0120E+01  7.2181E+01  4.2159E+01  5.7919E+01
 EBVSHRINKSD(%)  6.2023E-01  2.8935E+01  4.8658E+01  2.5777E+01  3.3619E+01
 EBVSHRINKVR(%)  1.2366E+00  4.9497E+01  7.3640E+01  4.4910E+01  5.5935E+01
 RELATIVEINF(%)  9.8669E+01  3.9586E+00  5.1693E+00  4.5631E+00  1.1737E+01
 EPSSHRINKSD(%)  3.3177E+01
 EPSSHRINKVR(%)  5.5347E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1981.5306915840047     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1062.5921583793320     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1981.531       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  1.68E+00  6.22E-01  6.10E-01  1.14E+00  9.77E-01  6.38E-01  2.08E+00  1.18E+00  9.32E-01  1.37E+00
 


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
 #CPUT: Total CPU Time in Seconds,       74.822
Stop Time:
Sat Oct 23 21:59:58 CDT 2021
