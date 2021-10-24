Sun Oct 24 03:56:42 CDT 2021
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
$DATA ../../../../data/SD4/TD1/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1685.28915702581        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2612E+02 -3.4592E+01 -6.0956E+01  5.2478E+01  7.9258E+01  7.2617E+01 -2.5358E+00  1.1339E+01  2.3431E+01  1.4244E+01
            -1.7996E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1693.21592098318        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  9.9209E-01  1.1099E+00  1.1716E+00  9.4677E-01  1.0391E+00  8.5941E-01  1.0218E+00  9.3663E-01  9.0933E-01  8.9253E-01
             1.1060E+00
 PARAMETER:  9.2063E-02  2.0423E-01  2.5835E-01  4.5305E-02  1.3834E-01 -5.1511E-02  1.2153E-01  3.4538E-02  4.9494E-03 -1.3691E-02
             2.0079E-01
 GRADIENT:  -2.3551E+01 -1.3361E+01  1.0658E+01 -2.9693E+01 -1.4966E+01 -2.8033E+01 -5.0671E+00  1.9573E-01 -4.6628E+00 -5.0836E+00
             1.1619E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1693.92757908624        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.8773E-01  1.1419E+00  1.1669E+00  9.3611E-01  1.0615E+00  8.9661E-01  1.0557E+00  7.6325E-01  8.5933E-01  9.6539E-01
             1.1068E+00
 PARAMETER:  8.7653E-02  2.3271E-01  2.5437E-01  3.3982E-02  1.5969E-01 -9.1333E-03  1.5418E-01 -1.7016E-01 -5.1603E-02  6.4780E-02
             2.0147E-01
 GRADIENT:  -3.4183E+01 -8.9746E-01  1.1734E+01 -1.8309E+01 -1.0910E+01 -1.0137E+01 -3.3316E+00 -9.7161E-01 -1.1418E+01 -4.7165E-01
             1.0979E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1695.72345340288        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  1.0033E+00  1.1798E+00  9.0774E-01  9.1071E-01  9.7939E-01  9.2539E-01  1.0618E+00  4.2494E-01  9.3273E-01  9.0263E-01
             1.0585E+00
 PARAMETER:  1.0331E-01  2.6536E-01  3.2068E-03  6.4734E-03  7.9173E-02  2.2459E-02  1.5999E-01 -7.5581E-01  3.0356E-02 -2.4371E-03
             1.5681E-01
 GRADIENT:   4.9194E+00  3.1619E+00 -6.3708E-01  4.7348E+00 -2.4983E+00  1.2934E+00  6.5257E-01  3.2830E-01  9.8252E-01  1.5902E+00
            -7.7309E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1695.85374217022        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  1.0017E+00  1.3571E+00  7.9509E-01  7.9433E-01  1.0145E+00  9.2231E-01  9.5500E-01  2.7494E-01  1.0148E+00  8.9469E-01
             1.0647E+00
 PARAMETER:  1.0169E-01  4.0538E-01 -1.2930E-01 -1.3026E-01  1.1439E-01  1.9128E-02  5.3955E-02 -1.1912E+00  1.1473E-01 -1.1281E-02
             1.6274E-01
 GRADIENT:  -1.7888E+00  3.5611E+00  4.2129E-01  2.6006E+00 -2.3530E+00 -4.9167E-01 -4.3007E-01  1.4176E-01 -1.5224E-01  1.9794E-01
             8.6446E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1695.90165698857        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      888
 NPARAMETR:  1.0034E+00  1.3900E+00  7.7110E-01  7.6842E-01  1.0229E+00  9.2431E-01  9.4104E-01  1.1647E-01  1.0357E+00  8.9346E-01
             1.0625E+00
 PARAMETER:  1.0340E-01  4.2928E-01 -1.5994E-01 -1.6341E-01  1.2268E-01  2.1289E-02  3.9229E-02 -2.0501E+00  1.3509E-01 -1.2654E-02
             1.6063E-01
 GRADIENT:   2.4917E+00 -1.7899E+00  1.3365E+00 -2.8933E+00 -1.0087E+00  2.7474E-01 -6.3421E-02  2.0156E-02 -3.4285E-01 -6.5557E-01
            -3.5895E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1695.91728643117        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  1.0024E+00  1.3918E+00  7.6599E-01  7.6895E-01  1.0219E+00  9.2368E-01  9.4088E-01  2.5522E-02  1.0368E+00  8.9639E-01
             1.0627E+00
 PARAMETER:  1.0235E-01  4.3061E-01 -1.6659E-01 -1.6273E-01  1.2163E-01  2.0608E-02  3.9063E-02 -3.5682E+00  1.3614E-01 -9.3776E-03
             1.6078E-01
 GRADIENT:  -4.1168E-01 -1.5588E-01 -2.1494E-03  2.7041E-01  3.3701E-01 -2.8582E-02  3.0987E-02  1.1782E-03  3.3257E-02  1.2484E-02
             4.0981E-02

0ITERATION NO.:   34    OBJECTIVE VALUE:  -1695.91868057128        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1193
 NPARAMETR:  1.0034E+00  1.3917E+00  7.6438E-01  7.6830E-01  1.0212E+00  9.2418E-01  9.4059E-01  1.0000E-02  1.0368E+00  8.9546E-01
             1.0625E+00
 PARAMETER:  1.0343E-01  4.3055E-01 -1.6869E-01 -1.6358E-01  1.2098E-01  2.1149E-02  3.8757E-02 -4.9780E+00  1.3615E-01 -1.0419E-02
             1.6060E-01
 GRADIENT:   2.3977E+00 -1.0729E+00 -1.2929E-02 -4.7294E-01  3.1987E-01  1.7915E-01 -1.7372E-02  0.0000E+00  3.2286E-02  1.9299E-02
            -9.0767E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1193
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.0524E-04 -1.3814E-02 -2.9522E-04  7.9887E-03 -2.4543E-02
 SE:             2.9787E-02  2.2919E-02  1.2516E-04  2.2332E-02  2.2589E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8379E-01  5.4669E-01  1.8333E-02  7.2055E-01  2.7727E-01

 ETASHRINKSD(%)  2.0934E-01  2.3219E+01  9.9581E+01  2.5185E+01  2.4323E+01
 ETASHRINKVR(%)  4.1825E-01  4.1047E+01  9.9998E+01  4.4027E+01  4.2730E+01
 EBVSHRINKSD(%)  5.5979E-01  2.2718E+01  9.9620E+01  2.6465E+01  2.3116E+01
 EBVSHRINKVR(%)  1.1164E+00  4.0275E+01  9.9999E+01  4.5926E+01  4.0889E+01
 RELATIVEINF(%)  9.8579E+01  2.5658E+00  1.2892E-04  2.3505E+00  7.7895E+00
 EPSSHRINKSD(%)  4.2619E+01
 EPSSHRINKVR(%)  6.7074E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1695.9186805712761     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -960.76785400753795     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1695.919       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.39E+00  7.64E-01  7.68E-01  1.02E+00  9.24E-01  9.41E-01  1.00E-02  1.04E+00  8.95E-01  1.06E+00
 


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
 #CPUT: Total CPU Time in Seconds,       34.025
Stop Time:
Sun Oct 24 03:56:50 CDT 2021
