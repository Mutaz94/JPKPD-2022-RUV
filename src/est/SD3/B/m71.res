Sat Oct 23 21:20:33 CDT 2021
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
$DATA ../../../../data/SD3/B/dat71.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2035.95601129109        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9272E+02  3.5366E+00  6.8452E+01 -4.1094E+01 -2.1395E+01  3.9510E+01 -1.4647E+01 -1.1448E+01 -5.9222E+00 -2.6838E+01
            -1.2009E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2046.37328486018        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.6465E-01  9.8980E-01  8.2618E-01  1.0983E+00  9.3985E-01  9.7419E-01  1.1129E+00  1.0332E+00  1.0328E+00  1.1672E+00
             1.0205E+00
 PARAMETER:  6.4006E-02  8.9750E-02 -9.0940E-02  1.9376E-01  3.7968E-02  7.3850E-02  2.0700E-01  1.3267E-01  1.3227E-01  2.5462E-01
             1.2031E-01
 GRADIENT:   3.4516E+01  1.2820E+01 -2.3583E+01  6.1316E+01  1.9846E+01 -3.5370E+00 -5.2457E+00  1.3411E+01  1.2483E+01  1.2391E+01
             1.9762E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2049.88957594151        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.6351E-01  9.6244E-01  6.3364E-01  1.1082E+00  7.8628E-01  9.8827E-01  1.4279E+00  7.1052E-01  9.2931E-01  1.0009E+00
             9.6938E-01
 PARAMETER:  6.2831E-02  6.1718E-02 -3.5627E-01  2.0270E-01 -1.4044E-01  8.8201E-02  4.5623E-01 -2.4175E-01  2.6692E-02  1.0088E-01
             6.8904E-02
 GRADIENT:   2.8207E+01  2.7567E+01 -4.3424E+01  9.3597E+01  3.0101E+01  1.2005E+00  1.1377E+01  9.1348E+00  3.5595E+00  1.5328E+01
            -8.9246E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2055.89702398315        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.4700E-01  9.0222E-01  6.2580E-01  1.0715E+00  7.6124E-01  9.8454E-01  1.3889E+00  4.0255E-01  9.0802E-01  9.1533E-01
             9.8447E-01
 PARAMETER:  4.5547E-02 -2.8976E-03 -3.6872E-01  1.6910E-01 -1.7281E-01  8.4421E-02  4.2853E-01 -8.0994E-01  3.5114E-03  1.1527E-02
             8.4348E-02
 GRADIENT:  -1.0406E+01 -1.4597E+01 -9.7263E+00 -1.1102E+01  1.6565E+01  2.9676E-01  7.1605E-01  1.3279E+00 -1.7622E+00  3.9659E-01
             1.0301E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2056.49057069689        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.5246E-01  9.4962E-01  6.0154E-01  1.0520E+00  7.5664E-01  9.8522E-01  1.3319E+00  1.6832E-01  9.2813E-01  9.2711E-01
             9.8498E-01
 PARAMETER:  5.1289E-02  4.8304E-02 -4.0826E-01  1.5068E-01 -1.7887E-01  8.5114E-02  3.8661E-01 -1.6819E+00  2.5415E-02  2.4319E-02
             8.4867E-02
 GRADIENT:   1.5571E+00 -1.8283E-02  1.2739E+00 -1.9291E+00 -5.6909E-01  3.0593E-01  2.3294E-01  4.7117E-02 -1.1136E+00 -2.1685E-02
             2.2799E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2056.51216116010        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  9.5192E-01  9.5553E-01  5.9407E-01  1.0484E+00  7.5441E-01  9.8463E-01  1.3207E+00  6.6021E-02  9.3555E-01  9.2679E-01
             9.8494E-01
 PARAMETER:  5.0728E-02  5.4510E-02 -4.2076E-01  1.4730E-01 -1.8182E-01  8.4512E-02  3.7812E-01 -2.6178E+00  3.3382E-02  2.3969E-02
             8.4829E-02
 GRADIENT:   4.2766E-02 -6.7819E-02 -2.1843E-01 -5.9040E-02  1.2370E-01 -1.4916E-02  3.7251E-02  3.0145E-03 -6.5162E-03  7.8137E-02
             6.3104E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2056.51303294768        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1058             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5250E-01  9.5630E-01  5.9335E-01  1.0478E+00  7.5412E-01  9.8499E-01  1.3203E+00  1.0000E-02  9.3594E-01  9.2671E-01
             9.8494E-01
 PARAMETER:  5.1331E-02  5.5317E-02 -4.2197E-01  1.4671E-01 -1.8220E-01  8.4873E-02  3.7784E-01 -4.5960E+00  3.3794E-02  2.3889E-02
             8.4828E-02
 GRADIENT:   3.7624E+02  3.0305E+01  1.6907E+01  1.1006E+02  2.1167E+01  4.1624E+01  1.6485E+01  0.0000E+00  5.5094E+00  9.3848E-01
             9.6863E-01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -2056.51303294768        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     1113
 NPARAMETR:  9.5250E-01  9.5630E-01  5.9335E-01  1.0478E+00  7.5412E-01  9.8499E-01  1.3203E+00  1.0000E-02  9.3594E-01  9.2671E-01
             9.8494E-01
 PARAMETER:  5.1331E-02  5.5317E-02 -4.2197E-01  1.4671E-01 -1.8220E-01  8.4873E-02  3.7784E-01 -4.5960E+00  3.3794E-02  2.3889E-02
             8.4828E-02
 GRADIENT:   1.3986E+00  8.9588E-02  2.0205E-01 -3.5866E-01 -3.4665E-01  1.2178E-01  6.2152E-02  0.0000E+00  9.2071E-03  1.5129E-03
            -3.1398E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1113
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.8901E-04 -1.2339E-03 -4.8157E-04 -5.4065E-04 -8.4388E-03
 SE:             2.9875E-02  2.2068E-02  2.1232E-04  2.5740E-02  2.3003E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8961E-01  9.5541E-01  2.3324E-02  9.8324E-01  7.1373E-01

 ETASHRINKSD(%)  1.0000E-10  2.6070E+01  9.9289E+01  1.3768E+01  2.2937E+01
 ETASHRINKVR(%)  1.0000E-10  4.5343E+01  9.9995E+01  2.5640E+01  4.0614E+01
 EBVSHRINKSD(%)  3.4523E-01  2.5601E+01  9.9385E+01  1.4038E+01  2.2751E+01
 EBVSHRINKVR(%)  6.8926E-01  4.4648E+01  9.9996E+01  2.6105E+01  4.0327E+01
 RELATIVEINF(%)  9.9138E+01  7.1981E+00  9.4610E-04  1.2371E+01  8.8045E+00
 EPSSHRINKSD(%)  3.4292E+01
 EPSSHRINKVR(%)  5.6824E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2056.5130329476801     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1137.5744997430074     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2056.513       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.52E-01  9.56E-01  5.93E-01  1.05E+00  7.54E-01  9.85E-01  1.32E+00  1.00E-02  9.36E-01  9.27E-01  9.85E-01
 


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
 #CPUT: Total CPU Time in Seconds,       90.356
Stop Time:
Sat Oct 23 21:20:49 CDT 2021
