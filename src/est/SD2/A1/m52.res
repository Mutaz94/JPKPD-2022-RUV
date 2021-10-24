Sat Oct 23 17:27:23 CDT 2021
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
$DATA ../../../../data/SD2/A1/dat52.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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
 RAW OUTPUT FILE (FILE): m52.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2617.72473592999        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1498E+02  2.5110E+01  6.0647E+01 -1.3449E+01  9.4785E+01  4.4551E+01 -6.5646E+01 -4.2211E+01 -1.9585E+01 -3.3236E+01
            -5.4687E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2693.18573030452        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1503E+00  8.6734E-01  7.6354E-01  1.1284E+00  7.7643E-01  1.2622E+00  1.2757E+00  9.8646E-01  1.0717E+00  9.4533E-01
             1.2979E+00
 PARAMETER:  2.4001E-01 -4.2319E-02 -1.6980E-01  2.2079E-01 -1.5305E-01  3.3285E-01  3.4350E-01  8.6367E-02  1.6928E-01  4.3784E-02
             3.6074E-01
 GRADIENT:   7.8390E+02  4.9774E+01 -4.1437E+00  1.3351E+02  3.9378E+01  1.1145E+02  4.3562E+00  1.2396E+01  2.4167E+01  9.9753E+00
            -2.7616E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2703.19593833993        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      200
 NPARAMETR:  1.1508E+00  5.9041E-01  6.4190E-01  1.2973E+00  5.8305E-01  1.1918E+00  1.4385E+00  7.9425E-01  1.0572E+00  7.5773E-01
             1.2877E+00
 PARAMETER:  2.4047E-01 -4.2693E-01 -3.4332E-01  3.6029E-01 -4.3948E-01  2.7550E-01  4.6363E-01 -1.3035E-01  1.5560E-01 -1.7743E-01
             3.5289E-01
 GRADIENT:   2.1442E+02 -1.0022E+01  3.3759E+01  8.5512E+01 -7.0011E+00  3.5904E+01 -6.0667E+00  1.3402E+01  1.6801E+01  8.9875E+00
            -3.9440E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2708.71639907927        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.1404E+00  6.2133E-01  6.0245E-01  1.2692E+00  5.9565E-01  1.2209E+00  1.6114E+00  3.3796E-01  1.0657E+00  7.5635E-01
             1.2683E+00
 PARAMETER:  2.3141E-01 -3.7590E-01 -4.0675E-01  3.3837E-01 -4.1811E-01  2.9955E-01  5.7711E-01 -9.8484E-01  1.6367E-01 -1.7926E-01
             3.3771E-01
 GRADIENT:   1.9162E+02 -3.1918E+00 -7.6211E+00  7.9392E+01  6.7495E+01  4.6643E+01  9.3194E+00 -5.0828E-01  1.9942E+01  4.7094E+00
            -7.4951E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2740.43404142928        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      557
 NPARAMETR:  9.9629E-01  5.1924E-01  4.5956E-01  1.2627E+00  4.6258E-01  1.0400E+00  1.5074E+00  1.5389E-01  9.7645E-01  6.4444E-01
             1.3361E+00
 PARAMETER:  9.6284E-02 -5.5539E-01 -6.7747E-01  3.3328E-01 -6.7094E-01  1.3919E-01  5.1036E-01 -1.7715E+00  7.6172E-02 -3.3937E-01
             3.8979E-01
 GRADIENT:  -1.8868E+01  1.2738E+01 -1.3840E+01  5.6353E+01  1.9481E+01  2.0849E+01 -2.0427E+00  2.4908E-01 -9.0860E+00 -1.0491E+00
             1.7241E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2743.48497314851        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      733            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0059E+00  4.5989E-01  4.0764E-01  1.2290E+00  4.1272E-01  9.8270E-01  1.5509E+00  3.8936E-02  1.0080E+00  6.4128E-01
             1.3066E+00
 PARAMETER:  1.0591E-01 -6.7676E-01 -7.9738E-01  3.0616E-01 -7.8499E-01  8.2553E-02  5.3885E-01 -3.1458E+00  1.0796E-01 -3.4430E-01
             3.6742E-01
 GRADIENT:   2.3815E+02  7.5434E+01  7.5952E+01  2.1616E+02  2.6119E+02  1.8863E+01  1.5333E+01  5.9286E-02  1.0392E+01  4.1574E+00
             4.5637E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2743.49308040773        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  1.0053E+00  4.5932E-01  4.0805E-01  1.2306E+00  4.1295E-01  9.8298E-01  1.5497E+00  1.8287E-02  1.0071E+00  6.4469E-01
             1.3070E+00
 PARAMETER:  1.0530E-01 -6.7801E-01 -7.9636E-01  3.0752E-01 -7.8442E-01  8.2829E-02  5.3808E-01 -3.9016E+00  1.0709E-01 -3.3898E-01
             3.6772E-01
 GRADIENT:   3.1299E-01 -1.9495E-02  6.8141E-01  6.6802E-02 -9.4779E-01 -2.9443E-02  6.2185E-02  4.6615E-03 -3.5533E-02  2.5190E-02
            -8.8212E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2743.49471876431        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1088
 NPARAMETR:  1.0055E+00  4.6002E-01  4.0823E-01  1.2297E+00  4.1299E-01  9.8316E-01  1.5499E+00  1.0000E-02  1.0073E+00  6.4478E-01
             1.3071E+00
 PARAMETER:  1.0550E-01 -6.7649E-01 -7.9592E-01  3.0676E-01 -7.8433E-01  8.3017E-02  5.3818E-01 -4.5490E+00  1.0727E-01 -3.3885E-01
             3.6781E-01
 GRADIENT:   7.7416E-01  7.9931E-01  1.8305E+00 -1.5069E+00 -2.8572E+00  4.9573E-02  1.6344E-01  2.4998E-05  3.9818E-02 -4.9858E-03
            -3.1484E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2743.49471876431        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1110
 NPARAMETR:  1.0055E+00  4.6002E-01  4.0823E-01  1.2297E+00  4.1299E-01  9.8316E-01  1.5499E+00  1.0000E-02  1.0073E+00  6.4478E-01
             1.3071E+00
 PARAMETER:  1.0550E-01 -6.7649E-01 -7.9592E-01  3.0676E-01 -7.8433E-01  8.3017E-02  5.3818E-01 -4.5490E+00  1.0727E-01 -3.3885E-01
             3.6781E-01
 GRADIENT:   7.7416E-01  7.9931E-01  1.8305E+00 -1.5069E+00 -2.8572E+00  4.9573E-02  1.6344E-01  2.4998E-05  3.9818E-02 -4.9858E-03
            -3.1484E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1110
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8423E-04  1.0495E-02 -3.4125E-04 -3.4308E-03  6.2263E-03
 SE:             2.9784E-02  2.2849E-02  2.9496E-04  2.9125E-02  2.3741E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8435E-01  6.4601E-01  2.4730E-01  9.0623E-01  7.9312E-01

 ETASHRINKSD(%)  2.1836E-01  2.3453E+01  9.9012E+01  2.4288E+00  2.0466E+01
 ETASHRINKVR(%)  4.3625E-01  4.1406E+01  9.9990E+01  4.7987E+00  3.6743E+01
 EBVSHRINKSD(%)  5.1815E-01  2.0770E+01  9.9016E+01  2.6783E+00  2.1509E+01
 EBVSHRINKVR(%)  1.0336E+00  3.7226E+01  9.9990E+01  5.2848E+00  3.8391E+01
 RELATIVEINF(%)  9.8956E+01  2.4527E+01  1.2510E-03  9.0723E+01  6.1175E+00
 EPSSHRINKSD(%)  2.5429E+01
 EPSSHRINKVR(%)  4.4392E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2743.4947187643074     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1456.9807722777657     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2743.495       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  4.60E-01  4.08E-01  1.23E+00  4.13E-01  9.83E-01  1.55E+00  1.00E-02  1.01E+00  6.45E-01  1.31E+00
 


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
 #CPUT: Total CPU Time in Seconds,       68.970
Stop Time:
Sat Oct 23 17:27:36 CDT 2021
