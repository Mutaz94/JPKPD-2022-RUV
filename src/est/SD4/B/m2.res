Sun Oct 24 01:32:39 CDT 2021
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
$DATA ../../../../data/SD4/B/dat2.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1658.80916365530        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4522E+02 -1.6170E+01 -1.5549E+01  1.7964E+00 -2.7518E+01  3.4158E+01  1.1710E+01  1.9249E+01  3.3043E+01  7.4158E+00
            -4.2769E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1670.74304471869        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.7015E-01  1.0034E+00  1.2763E+00  1.0679E+00  1.1680E+00  9.0191E-01  8.8666E-01  7.7122E-01  8.2161E-01  1.1091E+00
             1.0861E+00
 PARAMETER:  6.9701E-02  1.0335E-01  3.4394E-01  1.6570E-01  2.5526E-01 -3.2397E-03 -2.0297E-02 -1.5978E-01 -9.6492E-02  2.0358E-01
             1.8260E-01
 GRADIENT:   3.1517E+02  4.2089E+01 -6.4320E+00  1.2092E+02  2.0537E+01 -1.1031E+01 -1.6695E+00  3.1952E+00 -1.0398E+01 -9.0127E+00
            -5.8552E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1672.48791395593        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      206
 NPARAMETR:  9.6283E-01  1.0276E+00  1.2685E+00  1.0449E+00  1.1770E+00  9.2542E-01  7.3918E-01  5.5175E-01  9.6708E-01  1.2198E+00
             1.1434E+00
 PARAMETER:  6.2120E-02  1.2727E-01  3.3787E-01  1.4391E-01  2.6298E-01  2.2495E-02 -2.0222E-01 -4.9466E-01  6.6522E-02  2.9865E-01
             2.3398E-01
 GRADIENT:  -4.3608E+01  3.4037E+00 -5.5699E+00  1.5177E+01 -2.9467E+00 -3.1959E+01  1.8208E+00  4.6318E-01  9.1138E+00  3.1231E+00
             2.1068E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1674.88757484211        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      383
 NPARAMETR:  9.7830E-01  1.0377E+00  1.4187E+00  1.0331E+00  1.2497E+00  9.9687E-01  6.1174E-01  7.5554E-01  9.6912E-01  1.2820E+00
             1.0951E+00
 PARAMETER:  7.8056E-02  1.3705E-01  4.4974E-01  1.3254E-01  3.2288E-01  9.6866E-02 -3.9145E-01 -1.8033E-01  6.8633E-02  3.4843E-01
             1.9082E-01
 GRADIENT:  -7.4424E-02 -6.5468E-02  4.8770E-02  1.1842E+00 -2.2388E-01  3.4467E-01  3.2220E-02  5.0068E-02 -4.1750E-01  4.1269E-01
            -3.2506E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1674.95334953461        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      559
 NPARAMETR:  9.7828E-01  1.1777E+00  1.3375E+00  9.4135E-01  1.2832E+00  1.0017E+00  5.1947E-01  7.0857E-01  1.0594E+00  1.2940E+00
             1.0969E+00
 PARAMETER:  7.8043E-02  2.6360E-01  3.9081E-01  3.9556E-02  3.4939E-01  1.0174E-01 -5.5495E-01 -2.4450E-01  1.5766E-01  3.5772E-01
             1.9244E-01
 GRADIENT:  -2.3928E+00  6.3864E+00  1.4766E+00  5.4874E+00 -1.9134E+00  1.8941E+00 -1.2343E+00 -3.8014E-01 -2.9109E+00 -4.9779E-01
            -3.3323E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1675.18400738388        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  9.8148E-01  1.3423E+00  1.1628E+00  8.3751E-01  1.2979E+00  9.9392E-01  6.1973E-01  5.8814E-01  1.1270E+00  1.2712E+00
             1.0958E+00
 PARAMETER:  8.1304E-02  3.9440E-01  2.5079E-01 -7.7319E-02  3.6074E-01  9.3903E-02 -3.7847E-01 -4.3079E-01  2.1959E-01  3.3994E-01
             1.9148E-01
 GRADIENT:   2.1258E+00  9.1727E+00  9.3660E-01  9.3928E+00 -2.5052E+00 -1.8298E+00  6.8127E-02 -4.1801E-02  8.9273E-02  2.4009E-01
            -6.3594E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1675.29345900685        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      914
 NPARAMETR:  9.8261E-01  1.5328E+00  1.0072E+00  7.1458E-01  1.3428E+00  9.9414E-01  6.0767E-01  5.1739E-01  1.2614E+00  1.2691E+00
             1.0963E+00
 PARAMETER:  8.2452E-02  5.2708E-01  1.0713E-01 -2.3606E-01  3.9476E-01  9.4123E-02 -3.9813E-01 -5.5897E-01  3.3226E-01  3.3830E-01
             1.9193E-01
 GRADIENT:   2.3316E+00  1.5650E+01 -5.1965E-01  1.3189E+01 -9.6662E-01 -2.2223E+00 -1.9743E-01  2.4117E-01  1.1502E-01 -4.8984E-01
            -9.2957E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1675.41520432948        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  9.8133E-01  1.6942E+00  8.8819E-01  6.0264E-01  1.3930E+00  1.0043E+00  5.9167E-01  4.3534E-01  1.4115E+00  1.2867E+00
             1.0997E+00
 PARAMETER:  8.1150E-02  6.2723E-01 -1.8573E-02 -4.0643E-01  4.3148E-01  1.0425E-01 -4.2480E-01 -7.3163E-01  4.4468E-01  3.5211E-01
             1.9501E-01
 GRADIENT:  -1.7651E+00  1.1349E+01  5.3888E-01  6.5626E+00 -2.0704E+00  1.5812E+00 -3.5993E-01  2.1391E-01 -7.7833E-01 -1.3768E-01
            -3.4043E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1675.50227394457        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1275
 NPARAMETR:  9.8167E-01  1.7250E+00  8.4717E-01  5.7589E-01  1.4041E+00  9.9957E-01  5.8810E-01  2.6838E-01  1.4602E+00  1.2925E+00
             1.1006E+00
 PARAMETER:  8.1499E-02  6.4523E-01 -6.5858E-02 -4.5184E-01  4.3941E-01  9.9569E-02 -4.3085E-01 -1.2153E+00  4.7857E-01  3.5660E-01
             1.9589E-01
 GRADIENT:  -1.2295E+00 -3.6139E+00 -5.9166E-01  1.5054E+00  1.9422E+00 -3.0914E-01 -2.2172E-01  3.6602E-02  3.0066E-01  4.9809E-01
             3.8075E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1675.52379900813        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1456             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8307E-01  1.7468E+00  8.1464E-01  5.5958E-01  1.3993E+00  1.0008E+00  5.9375E-01  5.9391E-02  1.4764E+00  1.2737E+00
             1.0996E+00
 PARAMETER:  8.2921E-02  6.5781E-01 -1.0501E-01 -4.8057E-01  4.3597E-01  1.0076E-01 -4.2129E-01 -2.7236E+00  4.8959E-01  3.4195E-01
             1.9490E-01
 GRADIENT:   3.3935E+02  5.0984E+02  9.4644E-01  7.5071E+01  1.6242E+01  3.2388E+01  1.0687E+01  3.5285E-03  1.4571E+01  1.9958E+00
             1.0369E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1675.52706792772        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1635
 NPARAMETR:  9.8279E-01  1.7474E+00  8.1166E-01  5.5955E-01  1.3989E+00  1.0007E+00  5.9412E-01  2.4055E-02  1.4772E+00  1.2789E+00
             1.1003E+00
 PARAMETER:  8.2637E-02  6.5812E-01 -1.0867E-01 -4.8062E-01  4.3571E-01  1.0073E-01 -4.2068E-01 -3.6274E+00  4.9013E-01  3.4599E-01
             1.9557E-01
 GRADIENT:   1.0285E+00 -5.1023E+00  1.9315E-03 -6.1433E-02  3.7875E-02  6.7103E-02  6.1133E-02  2.4902E-04  1.0781E-01  8.7773E-02
             7.4290E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1675.52729888602        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1809
 NPARAMETR:  9.8294E-01  1.7438E+00  8.1056E-01  5.5962E-01  1.3957E+00  1.0009E+00  5.9465E-01  1.0000E-02  1.4708E+00  1.2748E+00
             1.0998E+00
 PARAMETER:  8.2735E-02  6.5797E-01 -1.0894E-01 -4.8045E-01  4.3541E-01  1.0077E-01 -4.2061E-01 -4.5507E+00  4.8990E-01  3.4551E-01
             1.9543E-01
 GRADIENT:  -7.6259E-03  3.1199E-01  1.1213E-02  1.9105E-03  1.0329E-01 -2.5343E-03 -4.0976E-03  0.0000E+00  3.7889E-02  2.9160E-02
             8.7900E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1809
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.5058E-04 -3.8360E-02 -1.9466E-04  2.0141E-02 -4.3121E-02
 SE:             2.9806E-02  1.9032E-02  8.8546E-05  2.3416E-02  2.3152E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9597E-01  4.3845E-02  2.7924E-02  3.8972E-01  6.2535E-02

 ETASHRINKSD(%)  1.4447E-01  3.6242E+01  9.9703E+01  2.1552E+01  2.2437E+01
 ETASHRINKVR(%)  2.8873E-01  5.9349E+01  9.9999E+01  3.8460E+01  3.9840E+01
 EBVSHRINKSD(%)  5.0270E-01  3.4320E+01  9.9737E+01  2.3482E+01  1.9622E+01
 EBVSHRINKVR(%)  1.0029E+00  5.6862E+01  9.9999E+01  4.1449E+01  3.5394E+01
 RELATIVEINF(%)  9.8875E+01  2.6146E+00  1.6456E-04  4.0469E+00  2.1493E+01
 EPSSHRINKSD(%)  4.1174E+01
 EPSSHRINKVR(%)  6.5395E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1675.5272988860179     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -940.37647232227971     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1675.527       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.75E+00  8.11E-01  5.60E-01  1.40E+00  1.00E+00  5.94E-01  1.00E-02  1.48E+00  1.28E+00  1.10E+00
 


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
 #CPUT: Total CPU Time in Seconds,       54.357
Stop Time:
Sun Oct 24 01:32:50 CDT 2021
