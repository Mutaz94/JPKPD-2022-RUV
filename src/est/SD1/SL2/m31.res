Sat Oct 23 15:05:17 CDT 2021
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
$DATA ../../../../data/SD1/SL2/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1211.66590323069        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9516E+02  2.9201E+01  8.7282E+01  1.9000E+02  1.8175E+02  4.1081E+01 -7.2566E+01 -1.2461E+02 -5.9676E+01 -3.0411E+01
            -4.9908E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2835.12402566680        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0344E+00  1.3054E+00  9.4477E-01  8.5735E-01  1.0553E+00  8.9868E-01  9.8449E-01  8.0850E-01  8.9330E-01  1.0737E+00
             2.1198E+00
 PARAMETER:  1.3379E-01  3.6655E-01  4.3181E-02 -5.3907E-02  1.5386E-01 -6.8306E-03  8.4368E-02 -1.1257E-01 -1.2830E-02  1.7113E-01
             8.5134E-01
 GRADIENT:   2.6786E+02  1.4944E+02 -7.4699E-01  3.0169E+01 -3.5081E+01 -3.1277E+01  1.1568E+01  4.4573E+00 -2.7931E+01 -2.0900E+01
            -3.3568E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2852.39934647500        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0226E+00  1.2417E+00  8.9499E-01  9.0943E-01  1.0312E+00  9.3238E-01  8.9901E-01  1.6672E-01  9.9152E-01  1.0894E+00
             2.1923E+00
 PARAMETER:  1.2232E-01  3.1651E-01 -1.0940E-02  5.0608E-03  1.3070E-01  2.9984E-02 -6.4614E-03 -1.6914E+00  9.1483E-02  1.8564E-01
             8.8497E-01
 GRADIENT:   2.0785E+02  1.1433E+02 -1.7349E+01  7.9498E+01  8.5772E+00 -1.2889E+01  3.2716E+00  2.8402E-01 -1.0603E+01 -9.2115E+00
            -2.4186E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2866.38497699216        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      247
 NPARAMETR:  9.6717E-01  1.1522E+00  9.3432E-01  9.1345E-01  9.9857E-01  9.4005E-01  8.8700E-01  6.7065E-02  1.0139E+00  1.0908E+00
             2.3554E+00
 PARAMETER:  6.6622E-02  2.4170E-01  3.2060E-02  9.4731E-03  9.8569E-02  3.8181E-02 -1.9911E-02 -2.6021E+00  1.1381E-01  1.8693E-01
             9.5671E-01
 GRADIENT:  -4.4047E+01 -2.4635E+01  7.8431E-01 -8.4895E+00 -7.8035E+00 -1.3243E+01 -5.7662E-01  4.6467E-02 -7.1688E-01 -5.0981E+00
            -6.0596E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2873.00392238085        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      424
 NPARAMETR:  9.8408E-01  1.5497E+00  9.6031E-01  6.8603E-01  1.2983E+00  9.7931E-01  6.4795E-01  1.6920E-01  1.2675E+00  1.3526E+00
             2.4102E+00
 PARAMETER:  8.3952E-02  5.3805E-01  5.9496E-02 -2.7684E-01  3.6109E-01  7.9088E-02 -3.3395E-01 -1.6767E+00  3.3706E-01  4.0203E-01
             9.7973E-01
 GRADIENT:  -4.2914E+00  1.4631E+00  1.4391E+00  3.8041E+00 -1.3897E+00  3.0458E+00 -5.1300E-01  9.2219E-02 -7.5973E-01  1.2896E+00
             7.0765E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2873.21142963484        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      599
 NPARAMETR:  9.8581E-01  1.6515E+00  9.1304E-01  6.1698E-01  1.3677E+00  9.7194E-01  6.2313E-01  1.4925E-01  1.3787E+00  1.3846E+00
             2.4090E+00
 PARAMETER:  8.5706E-02  6.0171E-01  9.0258E-03 -3.8291E-01  4.1312E-01  7.1544E-02 -3.7300E-01 -1.8021E+00  4.2111E-01  4.2538E-01
             9.7922E-01
 GRADIENT:  -3.6834E-01 -6.1900E-01 -6.5707E-02 -2.4746E-01 -2.2540E-02  1.9202E-01  9.3168E-02  6.1070E-02  1.5177E-01 -7.2028E-02
             5.5297E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2873.23979637516        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      777
 NPARAMETR:  9.8598E-01  1.6520E+00  9.1144E-01  6.1710E-01  1.3677E+00  9.7139E-01  6.2355E-01  3.7128E-02  1.3771E+00  1.3848E+00
             2.4087E+00
 PARAMETER:  8.5879E-02  6.0198E-01  7.2689E-03 -3.8272E-01  4.1312E-01  7.0970E-02 -3.7233E-01 -3.1934E+00  4.1999E-01  4.2556E-01
             9.7907E-01
 GRADIENT:   1.0518E-02  2.5293E-01 -8.9866E-02  1.3760E-01  1.1208E-01 -2.5389E-02  6.4080E-02  3.7596E-03 -4.0354E-02 -3.5255E-02
             3.8364E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2873.24295145863        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      955
 NPARAMETR:  9.8593E-01  1.6473E+00  9.1358E-01  6.2001E-01  1.3644E+00  9.7141E-01  6.2319E-01  1.0313E-02  1.3735E+00  1.3835E+00
             2.4087E+00
 PARAMETER:  8.5829E-02  5.9912E-01  9.6107E-03 -3.7802E-01  4.1073E-01  7.0995E-02 -3.7290E-01 -4.4744E+00  4.1735E-01  4.2462E-01
             9.7907E-01
 GRADIENT:  -9.1829E-02 -6.5427E-02 -2.1723E-02  8.1856E-02  2.0900E-02 -7.9633E-03 -1.6867E-02  2.8609E-04 -3.6296E-03 -1.5401E-02
             2.0925E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2873.24295692918        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  9.8594E-01  1.6473E+00  9.1365E-01  6.1997E-01  1.3645E+00  9.7142E-01  6.2333E-01  1.0000E-02  1.3734E+00  1.3835E+00
             2.4087E+00
 PARAMETER:  8.5836E-02  5.9915E-01  9.6902E-03 -3.7808E-01  4.1075E-01  7.1000E-02 -3.7269E-01 -4.5648E+00  4.1730E-01  4.2462E-01
             9.7907E-01
 GRADIENT:  -7.5793E-02 -5.0069E-02 -6.5677E-03  6.5774E-02 -2.6145E-03 -5.9505E-03 -8.0908E-03  0.0000E+00 -7.5860E-03 -1.2806E-02
             1.8245E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1012
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6697E-03 -3.1568E-02 -5.4501E-05  1.8313E-02 -1.6957E-02
 SE:             2.9442E-02  1.9071E-02  6.9936E-05  2.3615E-02  2.6654E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5478E-01  9.7868E-02  4.3580E-01  4.3807E-01  5.2465E-01

 ETASHRINKSD(%)  1.3649E+00  3.6109E+01  9.9766E+01  2.0886E+01  1.0705E+01
 ETASHRINKVR(%)  2.7112E+00  5.9179E+01  9.9999E+01  3.7410E+01  2.0263E+01
 EBVSHRINKSD(%)  1.4729E+00  3.5990E+01  9.9753E+01  2.1723E+01  1.0025E+01
 EBVSHRINKVR(%)  2.9242E+00  5.9027E+01  9.9999E+01  3.8727E+01  1.9046E+01
 RELATIVEINF(%)  9.7025E+01  5.0259E+00  3.9898E-04  8.7775E+00  2.2946E+01
 EPSSHRINKSD(%)  1.6394E+01
 EPSSHRINKVR(%)  3.0100E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          898
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1650.4136056355921     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2873.2429569291767     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1222.8293512935845     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2873.243       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.65E+00  9.14E-01  6.20E-01  1.36E+00  9.71E-01  6.23E-01  1.00E-02  1.37E+00  1.38E+00  2.41E+00
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       66.136
Stop Time:
Sat Oct 23 15:05:30 CDT 2021
