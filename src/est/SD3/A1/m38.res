Sat Oct 23 21:46:22 CDT 2021
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
$DATA ../../../../data/SD3/A1/dat38.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1341.40789671843        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7108E+02  8.6673E+01  5.7690E+01  7.9723E+01  1.0472E+02  1.8356E+01 -9.5693E+01 -1.8638E+01 -7.2480E+01 -8.5356E+01
            -1.1974E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1678.45371189111        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1413E+00  9.2221E-01  9.8425E-01  1.0737E+00  8.8942E-01  1.3621E+00  1.2862E+00  9.0234E-01  1.3501E+00  1.1334E+00
             2.4502E+00
 PARAMETER:  2.3214E-01  1.9020E-02  8.4126E-02  1.7107E-01 -1.7184E-02  4.0901E-01  3.5169E-01 -2.7691E-03  4.0020E-01  2.2522E-01
             9.9617E-01
 GRADIENT:   3.3248E+02  1.8516E+01  2.8214E+00  2.5934E+01 -2.3343E+01  8.5966E+01  2.6464E+00  1.1774E+01  3.6228E+01  9.1127E+00
             1.4761E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1711.97133453477        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0676E+00  4.9971E-01  6.1252E-01  1.2652E+00  5.7998E-01  1.2931E+00  1.4006E+00  4.1149E-02  1.0266E+00  8.8689E-01
             2.3464E+00
 PARAMETER:  1.6541E-01 -5.9373E-01 -3.9018E-01  3.3524E-01 -4.4477E-01  3.5704E-01  4.3691E-01 -3.0906E+00  1.2621E-01 -2.0029E-02
             9.5289E-01
 GRADIENT:   2.1563E+02  8.3011E+00 -5.5538E+01  9.2765E+01  9.9857E+01  8.5145E+01 -1.0998E+01  6.3712E-02 -1.4408E+01  1.1901E+01
             1.3940E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1733.87209265063        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      294
 NPARAMETR:  9.9934E-01  4.5583E-01  6.1798E-01  1.2449E+00  5.6030E-01  1.1549E+00  1.5260E+00  4.2116E-02  1.0242E+00  8.0066E-01
             2.1368E+00
 PARAMETER:  9.9340E-02 -6.8563E-01 -3.8130E-01  3.1903E-01 -4.7929E-01  2.4398E-01  5.2268E-01 -3.0673E+00  1.2389E-01 -1.2231E-01
             8.5930E-01
 GRADIENT:   1.2010E+01 -6.8730E+00 -2.7080E+01 -4.0443E+01  5.7592E+01  3.2103E+01 -1.4972E+01  5.6314E-02 -1.8845E+01 -5.8716E+00
             6.9364E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1750.47266014560        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      471
 NPARAMETR:  9.8928E-01  2.4598E-01  5.5542E-01  1.3576E+00  4.5969E-01  1.0413E+00  2.7162E+00  4.6513E-02  1.0249E+00  7.9831E-01
             1.8852E+00
 PARAMETER:  8.9221E-02 -1.3025E+00 -4.8802E-01  4.0572E-01 -6.7721E-01  1.4047E-01  1.0992E+00 -2.9680E+00  1.2458E-01 -1.2526E-01
             7.3403E-01
 GRADIENT:   1.0342E+01  9.6670E+00  1.7879E+01  1.1742E+01 -2.7415E+01 -4.1401E+00 -1.3159E+00  6.6178E-02 -5.5035E+00 -2.9595E+00
            -2.0933E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1752.17657321734        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  9.8453E-01  1.7200E-01  5.0668E-01  1.3483E+00  4.2338E-01  1.0472E+00  3.5006E+00  2.9095E-02  1.0393E+00  7.7539E-01
             1.9302E+00
 PARAMETER:  8.4404E-02 -1.6603E+00 -5.7987E-01  3.9884E-01 -7.5947E-01  1.4615E-01  1.3529E+00 -3.4372E+00  1.3850E-01 -1.5439E-01
             7.5761E-01
 GRADIENT:   5.1956E+00  2.8512E+00  1.0122E+01 -1.5559E+01 -1.8254E+01 -8.1229E-01  7.5198E-01  2.8036E-02  2.9079E+00  2.4695E+00
             4.4696E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1752.89740663869        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      824
 NPARAMETR:  9.7963E-01  1.2423E-01  5.1734E-01  1.3821E+00  4.2897E-01  1.0433E+00  4.3514E+00  2.0684E-02  1.0092E+00  7.5457E-01
             1.9251E+00
 PARAMETER:  7.9423E-02 -1.9856E+00 -5.5906E-01  4.2359E-01 -7.4636E-01  1.4243E-01  1.5705E+00 -3.7784E+00  1.0914E-01 -1.8161E-01
             7.5498E-01
 GRADIENT:   1.3805E+00  7.8704E+00 -1.3707E+01 -9.7087E+00  1.5689E+01 -1.5927E+00  1.2750E+01  1.1902E-02 -4.9295E+00 -5.9273E+00
            -1.8551E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1752.93056944141        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  9.8052E-01  1.2656E-01  5.2263E-01  1.3837E+00  4.3168E-01  1.0489E+00  4.2976E+00  1.0000E-02  1.0108E+00  7.5997E-01
             1.9223E+00
 PARAMETER:  8.0328E-02 -1.9671E+00 -5.4888E-01  4.2476E-01 -7.4007E-01  1.4779E-01  1.5581E+00 -4.6723E+00  1.1075E-01 -1.7448E-01
             7.5352E-01
 GRADIENT:   1.6962E+00  6.5976E-01  7.0999E-01  9.3916E-01 -2.3421E+00  3.1051E-01 -2.8170E-02  0.0000E+00 -3.0829E-01 -1.4772E-01
            -3.8470E-01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1752.93056944141        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  9.8052E-01  1.2656E-01  5.2263E-01  1.3837E+00  4.3168E-01  1.0489E+00  4.2976E+00  1.0000E-02  1.0108E+00  7.5997E-01
             1.9223E+00
 PARAMETER:  8.0328E-02 -1.9671E+00 -5.4888E-01  4.2476E-01 -7.4007E-01  1.4779E-01  1.5581E+00 -4.6723E+00  1.1075E-01 -1.7448E-01
             7.5352E-01
 GRADIENT:   1.6962E+00  6.5976E-01  7.0999E-01  9.3916E-01 -2.3421E+00  3.1051E-01 -2.8170E-02  0.0000E+00 -3.0829E-01 -1.4772E-01
            -3.8470E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1024
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.8726E-04  3.2608E-02 -1.8985E-04 -1.8674E-02  7.1093E-03
 SE:             2.9541E-02  1.4471E-02  1.8919E-04  2.7815E-02  2.1653E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7334E-01  2.4240E-02  3.1561E-01  5.0200E-01  7.4267E-01

 ETASHRINKSD(%)  1.0334E+00  5.1520E+01  9.9366E+01  6.8155E+00  2.7459E+01
 ETASHRINKVR(%)  2.0561E+00  7.6497E+01  9.9996E+01  1.3167E+01  4.7378E+01
 EBVSHRINKSD(%)  1.1276E+00  6.2884E+01  9.9222E+01  5.4335E+00  2.2968E+01
 EBVSHRINKVR(%)  2.2424E+00  8.6224E+01  9.9994E+01  1.0572E+01  4.0660E+01
 RELATIVEINF(%)  9.6987E+01  6.7558E+00  4.6524E-04  4.6584E+01  4.5254E+00
 EPSSHRINKSD(%)  2.9867E+01
 EPSSHRINKVR(%)  5.0814E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1752.9305694414109     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -833.99203623673816     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1752.931       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.27E-01  5.23E-01  1.38E+00  4.32E-01  1.05E+00  4.30E+00  1.00E-02  1.01E+00  7.60E-01  1.92E+00
 


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
 #CPUT: Total CPU Time in Seconds,       86.300
Stop Time:
Sat Oct 23 21:46:36 CDT 2021
