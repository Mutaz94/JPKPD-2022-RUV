Sat Oct 23 19:21:16 CDT 2021
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
$DATA ../../../../data/SD2/SL3/dat47.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      799
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

 TOT. NO. OF OBS RECS:      699
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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1864.73333722540        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3025E+02 -1.1321E+02 -2.5514E+01  5.9658E+01  1.7134E+02  5.6338E+01 -4.6610E+01 -6.0030E+01 -9.1341E+01 -5.1309E+01
            -1.9583E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2409.02088090863        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.5998E-01  1.2981E+00  8.5473E-01  8.4920E-01  9.3388E-01  9.4581E-01  9.8313E-01  1.1441E+00  1.3888E+00  1.3686E+00
             1.8249E+00
 PARAMETER:  5.9158E-02  3.6088E-01 -5.6969E-02 -6.3455E-02  3.1593E-02  4.4284E-02  8.2985E-02  2.3465E-01  4.2844E-01  4.1381E-01
             7.0155E-01
 GRADIENT:   3.0145E+01  6.2782E+01  7.0427E+00 -1.9495E+01 -7.0964E+01  6.3134E+00  1.6579E+01  3.7163E+00  2.4150E+01  2.6027E+01
            -6.4392E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2422.98589465957        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      220
 NPARAMETR:  9.7410E-01  1.5068E+00  9.9074E-01  7.7349E-01  1.1806E+00  9.4424E-01  9.2229E-01  1.2330E+00  1.4736E+00  1.3834E+00
             1.8406E+00
 PARAMETER:  7.3756E-02  5.1002E-01  9.0701E-02 -1.5684E-01  2.6601E-01  4.2623E-02  1.9106E-02  3.0947E-01  4.8769E-01  4.2456E-01
             7.1008E-01
 GRADIENT:  -6.0650E+01 -8.7462E+00 -8.3243E-02  1.2831E+01 -5.3504E+00 -7.1350E+00  1.1733E+01 -5.9800E+00  5.9262E+00  3.7741E-01
            -5.1334E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2430.82548303112        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      397
 NPARAMETR:  1.0027E+00  1.6197E+00  1.5595E+00  7.0876E-01  1.3933E+00  9.5588E-01  7.5612E-01  2.6097E+00  1.5251E+00  1.5177E+00
             1.8831E+00
 PARAMETER:  1.0269E-01  5.8227E-01  5.4438E-01 -2.4423E-01  4.3165E-01  5.4877E-02 -1.7955E-01  1.0592E+00  5.2207E-01  5.1723E-01
             7.3293E-01
 GRADIENT:   9.7146E+00 -4.3230E+00  3.7271E-02  2.8690E+00 -1.2521E+00 -5.4334E-01  2.4962E+00  2.1304E+00  4.3951E-01 -5.3489E+00
            -4.3911E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2432.06782953804        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      572
 NPARAMETR:  9.9951E-01  1.9573E+00  1.0409E+00  4.7833E-01  1.5094E+00  9.5516E-01  7.1415E-01  2.3499E+00  1.9940E+00  1.6623E+00
             1.8815E+00
 PARAMETER:  9.9506E-02  7.7154E-01  1.4008E-01 -6.3745E-01  5.1172E-01  5.4125E-02 -2.3667E-01  9.5438E-01  7.9014E-01  6.0823E-01
             7.3209E-01
 GRADIENT:   1.5275E+00  5.4308E+00  5.6735E-01  3.2395E+00 -2.6145E-03 -9.2017E-01  4.6333E-01  4.7634E-01  1.2006E+00 -1.0375E+00
            -2.5427E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2432.08893948063        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      747
 NPARAMETR:  9.9983E-01  2.0101E+00  9.1957E-01  4.4169E-01  1.5184E+00  9.5713E-01  7.1311E-01  2.2177E+00  2.0804E+00  1.6747E+00
             1.8832E+00
 PARAMETER:  9.9829E-02  7.9817E-01  1.6150E-02 -7.1715E-01  5.1764E-01  5.6180E-02 -2.3811E-01  8.9647E-01  8.3258E-01  6.1565E-01
             7.3300E-01
 GRADIENT:   2.0432E+00  6.1825E+00  1.7360E-01  2.7788E+00 -2.6068E-01 -1.5340E-01 -1.9784E-01  1.8653E-01  4.3374E-01 -6.2103E-01
            -1.5464E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2432.11393166288        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      916
 NPARAMETR:  9.9911E-01  2.0055E+00  9.1381E-01  4.3939E-01  1.5180E+00  9.5749E-01  7.1322E-01  2.1981E+00  2.0810E+00  1.6792E+00
             1.8846E+00
 PARAMETER:  9.9107E-02  7.9588E-01  9.8725E-03 -7.2237E-01  5.1737E-01  5.6562E-02 -2.3797E-01  8.8758E-01  8.3287E-01  6.1834E-01
             7.3370E-01
 GRADIENT:   4.1725E-01 -2.6986E+00  1.2299E-02  2.4855E-01  2.8309E-01  2.9467E-02  1.3928E-02  9.0238E-02  2.7361E-01  2.6281E-01
            -4.1771E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      916
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2175E-03 -5.4633E-02 -2.4097E-02  4.1451E-02 -3.8889E-02
 SE:             2.9640E-02  2.0818E-02  1.1152E-02  2.3078E-02  2.4663E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6723E-01  8.6832E-03  3.0710E-02  7.2475E-02  1.1484E-01

 ETASHRINKSD(%)  7.0294E-01  3.0256E+01  6.2640E+01  2.2686E+01  1.7374E+01
 ETASHRINKVR(%)  1.4009E+00  5.1358E+01  8.6042E+01  4.0225E+01  3.1730E+01
 EBVSHRINKSD(%)  1.0200E+00  2.8743E+01  6.6063E+01  2.6941E+01  1.2197E+01
 EBVSHRINKVR(%)  2.0297E+00  4.9225E+01  8.8483E+01  4.6625E+01  2.2907E+01
 RELATIVEINF(%)  9.7945E+01  9.0547E+00  6.3494E+00  1.0049E+01  4.3484E+01
 EPSSHRINKSD(%)  2.3286E+01
 EPSSHRINKVR(%)  4.1150E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          699
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1284.6760694201323     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2432.1139316628814     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1147.4378622427491     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2432.114       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  2.01E+00  9.14E-01  4.39E-01  1.52E+00  9.57E-01  7.13E-01  2.20E+00  2.08E+00  1.68E+00  1.88E+00
 


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
 #CPUT: Total CPU Time in Seconds,       83.797
Stop Time:
Sat Oct 23 19:21:29 CDT 2021
