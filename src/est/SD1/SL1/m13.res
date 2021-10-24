Sat Oct 23 14:44:27 CDT 2021
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
$DATA ../../../../data/SD1/SL1/dat13.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m13.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3114.38228198190        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4342E+02  8.9948E+01 -2.9151E+01  4.2130E+01  7.7094E+01  3.8576E+01 -4.4611E+01 -6.5100E+01 -1.3359E+01  1.6776E+01
            -1.3775E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3382.13045700552        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.6274E-01  9.4774E-01  1.1690E+00  1.0740E+00  9.9802E-01  9.9570E-01  1.2540E+00  1.1128E+00  8.9491E-01  7.5827E-01
             1.5309E+00
 PARAMETER:  6.2023E-02  4.6326E-02  2.5615E-01  1.7139E-01  9.8018E-02  9.5688E-02  3.2634E-01  2.0691E-01 -1.1032E-02 -1.7672E-01
             5.2585E-01
 GRADIENT:   1.0295E+02  1.7617E+01  2.1943E+01  1.0409E+02  2.6971E+01  1.3068E+01  3.2149E+00 -1.2520E+01 -1.1199E+01  2.8661E-02
            -3.5475E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3383.81423798814        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      201
 NPARAMETR:  9.7040E-01  9.3353E-01  1.1238E+00  1.0494E+00  9.8706E-01  9.9485E-01  1.4975E+00  1.2802E+00  1.0196E+00  6.1854E-01
             1.5157E+00
 PARAMETER:  6.9955E-02  3.1214E-02  2.1672E-01  1.4826E-01  8.6975E-02  9.4833E-02  5.0382E-01  3.4701E-01  1.1944E-01 -3.8040E-01
             5.1589E-01
 GRADIENT:  -6.3683E+01 -1.3217E+01 -8.9509E+00 -6.1464E+00  2.9766E+01 -7.5452E+00  1.2935E+01 -2.8135E+00  7.3314E+00 -3.0337E+00
            -4.5380E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3386.71798585247        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  9.9793E-01  1.0638E+00  1.4386E+00  1.0194E+00  1.1492E+00  1.0353E+00  1.2068E+00  1.8545E+00  9.2181E-01  8.6384E-01
             1.5223E+00
 PARAMETER:  9.7929E-02  1.6188E-01  4.6365E-01  1.1926E-01  2.3905E-01  1.3470E-01  2.8797E-01  7.1760E-01  1.8585E-02 -4.6363E-02
             5.2020E-01
 GRADIENT:  -2.0514E+00  4.3934E+00 -4.9868E+00  1.4241E+00  6.0788E+00  1.0008E+01 -1.5005E+00  5.0477E-01  8.5090E-01 -2.6722E-01
             2.6175E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3386.92974611682        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  9.9883E-01  1.0201E+00  1.3773E+00  1.0376E+00  1.0994E+00  1.0098E+00  1.2643E+00  1.7224E+00  9.2819E-01  8.0796E-01
             1.5199E+00
 PARAMETER:  9.8834E-02  1.1989E-01  4.2011E-01  1.3687E-01  1.9476E-01  1.0971E-01  3.3451E-01  6.4374E-01  2.5487E-02 -1.1324E-01
             5.1867E-01
 GRADIENT:  -1.3754E-01  2.3037E-01 -3.8360E-01  7.2051E-01 -3.3405E-01  5.3571E-01 -1.6626E-02  1.2032E-01  1.4557E-01  5.9045E-02
            -2.8465E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3387.36815528518        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  9.9868E-01  9.9041E-01  1.6449E+00  1.0671E+00  1.1415E+00  1.0091E+00  1.3011E+00  2.1469E+00  8.6292E-01  7.8834E-01
             1.5121E+00
 PARAMETER:  9.8683E-02  9.0367E-02  5.9765E-01  1.6492E-01  2.3234E-01  1.0906E-01  3.6322E-01  8.6401E-01 -4.7431E-02 -1.3782E-01
             5.1353E-01
 GRADIENT:  -1.0976E-01 -4.5183E-01  1.2268E+00  3.6154E-01  5.2893E-01  3.0077E-01  1.3300E-01 -1.5958E+00  2.0189E-02 -6.1689E-01
             8.6097E-01

0ITERATION NO.:   27    OBJECTIVE VALUE:  -3387.38473893087        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      791
 NPARAMETR:  9.9858E-01  9.9469E-01  1.6440E+00  1.0651E+00  1.1446E+00  1.0088E+00  1.2955E+00  2.1755E+00  8.6264E-01  7.9871E-01
             1.5113E+00
 PARAMETER:  9.8583E-02  9.4679E-02  5.9712E-01  1.6307E-01  2.3508E-01  1.0876E-01  3.5892E-01  8.7725E-01 -4.7763E-02 -1.2476E-01
             5.1299E-01
 GRADIENT:  -3.3208E-01  1.2865E-01  3.5597E-02  6.9518E-02 -1.0728E-01  1.7589E-01  1.7757E-01  8.7752E-02  1.1430E-01  1.0805E-02
             5.8101E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      791
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3518E-03 -1.1581E-02 -2.9957E-02  6.8213E-03 -3.8816E-02
 SE:             2.9784E-02  2.3505E-02  2.2307E-02  2.4739E-02  1.9728E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6380E-01  6.2222E-01  1.7930E-01  7.8276E-01  4.9125E-02

 ETASHRINKSD(%)  2.1925E-01  2.1255E+01  2.5269E+01  1.7120E+01  3.3907E+01
 ETASHRINKVR(%)  4.3802E-01  3.7992E+01  4.4153E+01  3.1309E+01  5.6318E+01
 EBVSHRINKSD(%)  5.6519E-01  2.1238E+01  2.6207E+01  1.8478E+01  3.3226E+01
 EBVSHRINKVR(%)  1.1272E+00  3.7965E+01  4.5546E+01  3.3542E+01  5.5413E+01
 RELATIVEINF(%)  9.8866E+01  2.9846E+01  3.6106E+01  3.6000E+01  2.1163E+01
 EPSSHRINKSD(%)  2.0185E+01
 EPSSHRINKVR(%)  3.6296E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          898
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1650.4136056355921     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3387.3847389308744     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1736.9711332952822     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3387.385       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  9.95E-01  1.64E+00  1.07E+00  1.14E+00  1.01E+00  1.30E+00  2.18E+00  8.63E-01  7.99E-01  1.51E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       52.449
Stop Time:
Sat Oct 23 14:44:36 CDT 2021
