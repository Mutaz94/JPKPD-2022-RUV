Sun Oct 24 02:52:03 CDT 2021
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
$DATA ../../../../data/SD4/S1/dat89.csv ignore=@
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
Current Date:       24 OCT 2021
Days until program expires : 175
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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1673.17554595834        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0642E+02 -1.8455E+01 -1.8221E+01  2.3379E+01  5.1734E+01  6.4018E+01  2.9641E+00  2.4490E+00  2.3480E+01 -8.3613E+00
             1.2967E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1677.67892244845        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9621E-01  1.0419E+00  1.0051E+00  1.0066E+00  9.7931E-01  8.6612E-01  9.9455E-01  9.9605E-01  9.2233E-01  1.0304E+00
             9.6037E-01
 PARAMETER:  9.6207E-02  1.4100E-01  1.0507E-01  1.0663E-01  7.9094E-02 -4.3727E-02  9.4536E-02  9.6040E-02  1.9146E-02  1.2992E-01
             5.9566E-02
 GRADIENT:   3.6121E+00  4.1707E+00  4.3341E+00  4.4220E-01 -7.0405E+00 -1.5040E+01 -1.2459E+00  7.7747E-02  1.9061E+00 -1.9038E+00
            -6.7254E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1678.12098707614        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.9337E-01  8.9001E-01  1.0443E+00  1.1134E+00  9.3253E-01  8.8572E-01  1.2322E+00  8.9513E-01  7.8293E-01  1.0186E+00
             9.7707E-01
 PARAMETER:  9.3350E-02 -1.6528E-02  1.4339E-01  2.0741E-01  3.0143E-02 -2.1349E-02  3.0882E-01 -1.0789E-02 -1.4472E-01  1.1847E-01
             7.6802E-02
 GRADIENT:  -2.3554E+00  1.9915E+01  7.9894E+00  1.8076E+01 -6.8214E+00 -5.6572E+00 -1.0303E+00 -2.4922E+00 -8.6899E+00  1.0337E+00
            -2.8272E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1679.49372489223        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.9429E-01  7.4613E-01  8.6720E-01  1.1931E+00  7.8707E-01  9.0906E-01  1.4525E+00  6.3902E-01  7.6573E-01  8.5475E-01
             9.7423E-01
 PARAMETER:  9.4276E-02 -1.9286E-01 -4.2480E-02  2.7652E-01 -1.3944E-01  4.6569E-03  4.7330E-01 -3.4783E-01 -1.6692E-01 -5.6951E-02
             7.3888E-02
 GRADIENT:  -1.4385E+00  2.0825E+01  3.7943E+00  4.3791E+01 -4.6208E-01  3.9543E+00 -2.3469E+00 -9.9063E-01  2.3149E-02 -3.3162E+00
            -8.7540E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1683.25409771415        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.9699E-01  5.1257E-01  5.5486E-01  1.2515E+00  5.2752E-01  8.9507E-01  1.9707E+00  2.3125E-01  6.9814E-01  6.0168E-01
             9.7790E-01
 PARAMETER:  9.6987E-02 -5.6832E-01 -4.8905E-01  3.2434E-01 -5.3956E-01 -1.0849E-02  7.7838E-01 -1.3643E+00 -2.5933E-01 -4.0802E-01
             7.7655E-02
 GRADIENT:  -9.1257E-03  1.4319E+01  9.1878E+00  3.0967E+01 -2.0531E+01 -2.9154E+00 -7.7453E-01  1.3674E-01 -3.5860E+00 -1.2987E+00
            -1.9526E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1683.86955931974        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  9.9548E-01  4.3144E-01  5.6636E-01  1.2796E+00  5.2058E-01  9.0029E-01  2.2112E+00  1.7294E-01  6.9652E-01  6.2801E-01
             9.8280E-01
 PARAMETER:  9.5471E-02 -7.4063E-01 -4.6852E-01  3.4654E-01 -5.5281E-01 -5.0342E-03  8.9354E-01 -1.6548E+00 -2.6166E-01 -3.6519E-01
             8.2655E-02
 GRADIENT:  -6.3757E-02  7.3714E-01  1.3085E+00  8.6223E-01 -1.0061E+00 -2.3615E-01 -7.2953E-02 -1.9856E-01  1.0705E-01 -3.8558E-01
            -6.4551E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1683.88569672530        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  9.9468E-01  3.9593E-01  5.7513E-01  1.2974E+00  5.1780E-01  9.0209E-01  2.3404E+00  2.3102E-01  6.9060E-01  6.3668E-01
             9.8229E-01
 PARAMETER:  9.4662E-02 -8.2651E-01 -4.5316E-01  3.6036E-01 -5.5816E-01 -3.0458E-03  9.5032E-01 -1.3653E+00 -2.7019E-01 -3.5149E-01
             8.2128E-02
 GRADIENT:   7.0223E-02 -8.3778E-01 -2.0992E-01 -4.2153E+00 -1.3126E+00  7.7040E-01 -1.2070E-01 -8.0799E-02  3.3619E-02  1.0302E+00
             2.1263E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1683.88856248390        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1230
 NPARAMETR:  9.9421E-01  3.8728E-01  5.8897E-01  1.3061E+00  5.2415E-01  9.0140E-01  2.3816E+00  2.7649E-01  6.8900E-01  6.4199E-01
             9.8150E-01
 PARAMETER:  9.4197E-02 -8.4861E-01 -4.2937E-01  3.6706E-01 -5.4598E-01 -3.8070E-03  9.6778E-01 -1.1856E+00 -2.7252E-01 -3.4318E-01
             8.1331E-02
 GRADIENT:   7.1896E-02 -6.1885E-01 -5.1619E-01 -3.2923E+00 -1.2073E+00  5.8336E-01 -1.9979E-02  5.0662E-02  6.4058E-02  1.0482E+00
             2.7337E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1683.89227752620        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1406
 NPARAMETR:  9.9446E-01  3.8788E-01  5.9641E-01  1.3084E+00  5.2847E-01  8.9975E-01  2.3974E+00  2.8828E-01  6.8875E-01  6.3813E-01
             9.8062E-01
 PARAMETER:  9.4447E-02 -8.4705E-01 -4.1683E-01  3.6878E-01 -5.3777E-01 -5.6421E-03  9.7437E-01 -1.1438E+00 -2.7288E-01 -3.4922E-01
             8.0428E-02
 GRADIENT:   1.2156E+00  3.5453E-01  2.6546E+00 -2.1112E+00 -3.1743E+00 -1.1031E-01  6.8649E-01 -1.5860E-01 -1.7354E-01 -5.8208E-01
            -3.4543E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1683.89666442343        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1584
 NPARAMETR:  9.9447E-01  3.8868E-01  5.9617E-01  1.3075E+00  5.2862E-01  9.0009E-01  2.3934E+00  2.9230E-01  6.8913E-01  6.3829E-01
             9.8076E-01
 PARAMETER:  9.4454E-02 -8.4500E-01 -4.1724E-01  3.6815E-01 -5.3748E-01 -5.2634E-03  9.7271E-01 -1.1300E+00 -2.7232E-01 -3.4896E-01
             8.0572E-02
 GRADIENT:   1.1776E+00  1.3523E-01  1.6407E+00 -3.1201E+00 -2.2368E+00  3.6576E-02  7.0574E-01 -6.2718E-02 -3.5065E-02 -2.3372E-01
            -1.3611E-01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1683.89764929323        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1682
 NPARAMETR:  9.9447E-01  3.8877E-01  5.9590E-01  1.3074E+00  5.2877E-01  9.0009E-01  2.3925E+00  2.9562E-01  6.8928E-01  6.3900E-01
             9.8093E-01
 PARAMETER:  9.4458E-02 -8.4477E-01 -4.1769E-01  3.6803E-01 -5.3720E-01 -5.2566E-03  9.7233E-01 -1.1187E+00 -2.7211E-01 -3.4785E-01
             8.0741E-02
 GRADIENT:  -3.9444E-02 -3.3481E-01 -3.0445E-01  9.1072E-01  3.6380E-01 -2.4152E-03  1.5787E-01  1.7957E-02 -1.7631E-03  1.5099E-01
             7.6597E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1682
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1246E-03  3.9907E-02 -1.7568E-02 -2.8088E-02  1.6316E-02
 SE:             2.9843E-02  2.0270E-02  8.1391E-03  2.5497E-02  2.1519E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6994E-01  4.8981E-02  3.0890E-02  2.7063E-01  4.4832E-01

 ETASHRINKSD(%)  2.2151E-02  3.2092E+01  7.2733E+01  1.4582E+01  2.7909E+01
 ETASHRINKVR(%)  4.4298E-02  5.3885E+01  9.2565E+01  2.7038E+01  4.8030E+01
 EBVSHRINKSD(%)  5.0956E-01  3.5025E+01  7.3357E+01  1.3025E+01  2.4921E+01
 EBVSHRINKVR(%)  1.0165E+00  5.7783E+01  9.2901E+01  2.4353E+01  4.3631E+01
 RELATIVEINF(%)  9.8272E+01  9.2338E+00  4.4940E-01  2.2340E+01  3.1758E+00
 EPSSHRINKSD(%)  4.4633E+01
 EPSSHRINKVR(%)  6.9345E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1683.8976492932270     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -948.74682272948883     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.69
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1683.898       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  3.89E-01  5.96E-01  1.31E+00  5.29E-01  9.00E-01  2.39E+00  2.96E-01  6.89E-01  6.39E-01  9.81E-01
 


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
 #CPUT: Total CPU Time in Seconds,       49.089
Stop Time:
Sun Oct 24 02:52:13 CDT 2021
