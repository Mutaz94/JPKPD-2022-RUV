Sat Oct 23 22:58:53 CDT 2021
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
$DATA ../../../../data/SD3/S1/dat31.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2077.00935930756        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8936E+02 -8.6516E+00 -4.0370E+01  6.4727E+01  9.0568E+01  3.0187E+01  1.0202E+01 -1.1122E+00  1.2741E+01  6.3013E+00
            -3.5984E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2082.34174066213        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7297E-01  1.0298E+00  1.0258E+00  1.0111E+00  9.6188E-01  1.0194E+00  9.6823E-01  1.0118E+00  9.8336E-01  9.5802E-01
             1.0537E+00
 PARAMETER:  7.2601E-02  1.2936E-01  1.2548E-01  1.1107E-01  6.1134E-02  1.1919E-01  6.7711E-02  1.1169E-01  8.3221E-02  5.7109E-02
             1.5233E-01
 GRADIENT:   4.4766E+00  7.5730E+00 -4.2668E+00  1.7194E+01  8.0404E+00  3.1793E-01  5.6186E+00 -3.5172E+00  9.1060E-01  6.0065E+00
             6.4160E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2083.54701146547        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  9.7400E-01  1.0916E+00  1.1412E+00  9.6905E-01  1.0171E+00  1.0251E+00  6.8563E-01  1.4069E+00  1.0432E+00  9.5278E-01
             1.0233E+00
 PARAMETER:  7.3658E-02  1.8766E-01  2.3209E-01  6.8564E-02  1.1695E-01  1.2480E-01 -2.7742E-01  4.4140E-01  1.4227E-01  5.1629E-02
             1.2299E-01
 GRADIENT:   6.5771E+00  8.7508E+00 -9.5798E-03  6.1574E+00 -2.9099E+00  2.2334E+00 -1.6094E+00  7.0582E-01 -9.5103E+00  7.9734E-02
            -1.5040E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2083.90947095715        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  9.7103E-01  1.1072E+00  1.1571E+00  9.5459E-01  1.0325E+00  1.0193E+00  6.4405E-01  1.4217E+00  1.1011E+00  9.6575E-01
             1.0423E+00
 PARAMETER:  7.0603E-02  2.0183E-01  2.4594E-01  5.3526E-02  1.3195E-01  1.1910E-01 -3.3998E-01  4.5185E-01  1.9635E-01  6.5154E-02
             1.4142E-01
 GRADIENT:  -1.6349E-01  6.3378E-01  9.1580E-02  1.0037E+00 -1.6999E-01  2.2418E-02 -1.0770E-01 -2.2366E-02 -2.8998E-01 -4.8489E-02
             3.3793E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2083.91906972903        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.7176E-01  1.1954E+00  1.0783E+00  8.9670E-01  1.0418E+00  1.0200E+00  6.7116E-01  1.4133E+00  1.1398E+00  9.6056E-01
             1.0418E+00
 PARAMETER:  7.1349E-02  2.7845E-01  1.7541E-01 -9.0288E-03  1.4096E-01  1.1976E-01 -2.9875E-01  4.4595E-01  2.3084E-01  5.9762E-02
             1.4098E-01
 GRADIENT:   2.5188E-03  6.5873E-01  4.4912E-01  2.5757E-01 -7.0958E-01  4.7610E-02 -5.4499E-02 -5.3209E-02 -4.7921E-02 -2.5339E-01
            -1.4526E-01

0ITERATION NO.:   23    OBJECTIVE VALUE:  -2083.92018075258        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      803
 NPARAMETR:  9.7220E-01  1.1949E+00  1.0771E+00  8.9660E-01  1.0421E+00  1.0201E+00  6.7255E-01  1.4141E+00  1.1400E+00  9.6188E-01
             1.0419E+00
 PARAMETER:  7.1808E-02  2.7804E-01  1.7426E-01 -9.1437E-03  1.4121E-01  1.1991E-01 -2.9668E-01  4.4646E-01  2.3103E-01  6.1130E-02
             1.4106E-01
 GRADIENT:   9.8791E-01 -5.5151E-01 -1.3799E-01 -1.8454E-01  2.1673E-01  1.0840E-01  4.9208E-02  1.3470E-01  1.3855E-01 -4.1883E-03
            -7.4946E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      803
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1197E-04 -3.2110E-02 -3.2460E-02  1.2657E-02 -3.8312E-02
 SE:             2.9866E-02  1.4870E-02  1.5475E-02  2.6515E-02  2.1901E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7831E-01  3.0821E-02  3.5945E-02  6.3311E-01  8.0231E-02

 ETASHRINKSD(%)  1.0000E-10  5.0184E+01  4.8156E+01  1.1170E+01  2.6630E+01
 ETASHRINKVR(%)  1.0000E-10  7.5184E+01  7.3122E+01  2.1093E+01  4.6168E+01
 EBVSHRINKSD(%)  3.6169E-01  5.0464E+01  5.0724E+01  1.1385E+01  2.4843E+01
 EBVSHRINKVR(%)  7.2207E-01  7.5462E+01  7.5719E+01  2.1474E+01  4.3514E+01
 RELATIVEINF(%)  9.8957E+01  1.2137E+00  3.8892E+00  4.6928E+00  1.4156E+01
 EPSSHRINKSD(%)  3.4276E+01
 EPSSHRINKVR(%)  5.6803E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2083.9201807525806     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1164.9816475479079     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2083.920       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.72E-01  1.19E+00  1.08E+00  8.97E-01  1.04E+00  1.02E+00  6.73E-01  1.41E+00  1.14E+00  9.62E-01  1.04E+00
 


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
 #CPUT: Total CPU Time in Seconds,       64.960
Stop Time:
Sat Oct 23 22:59:05 CDT 2021
