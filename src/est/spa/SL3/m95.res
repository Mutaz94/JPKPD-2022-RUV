Sat Sep 18 13:09:27 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat95.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1572.06001701983        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8284E+01 -8.4066E+01 -5.3234E+01 -9.2693E+01  1.1717E+02  1.2769E+01 -2.0485E+01  1.3851E+01 -8.3144E+01 -9.3837E+00
            -7.7914E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1593.04134977701        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0227E+00  9.0001E-01  1.0954E+00  1.0574E+00  9.5312E-01  9.5854E-01  1.0372E+00  8.5495E-01  1.3942E+00  9.1826E-01
             1.1869E+00
 PARAMETER:  1.2242E-01 -5.3482E-03  1.9109E-01  1.5583E-01  5.1991E-02  5.7653E-02  1.3653E-01 -5.6713E-02  4.3232E-01  1.4721E-02
             2.7133E-01
 GRADIENT:   8.8006E+01 -5.8732E+01 -2.3840E+01 -5.0113E+01  6.5480E+01 -4.7475E+00  4.2261E+00  7.8828E+00  2.5996E+01 -1.3093E+01
             7.5227E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1600.74733599436        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0161E+00  6.3138E-01  9.7957E-01  1.3101E+00  7.8889E-01  9.7052E-01  9.1581E-01  2.3380E-01  1.2065E+00  1.0453E+00
             1.2043E+00
 PARAMETER:  1.1593E-01 -3.5984E-01  7.9354E-02  3.7014E-01 -1.3713E-01  7.0077E-02  1.2053E-02 -1.3533E+00  2.8773E-01  1.4430E-01
             2.8590E-01
 GRADIENT:   6.8928E+01 -9.9827E-01 -3.2653E+01  4.9055E+01  2.6526E+01  1.8841E+00  1.4754E+00  7.8823E-01  2.3982E+01  7.6764E+00
             2.1765E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1603.21365311960        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.9368E-01  6.6950E-01  9.4964E-01  1.2535E+00  7.7607E-01  9.6538E-01  1.4940E+00  2.4685E-01  1.0993E+00  9.4448E-01
             1.1444E+00
 PARAMETER:  9.3664E-02 -3.0122E-01  4.8328E-02  3.2597E-01 -1.5351E-01  6.4761E-02  5.0146E-01 -1.2990E+00  1.9466E-01  4.2877E-02
             2.3489E-01
 GRADIENT:   1.5625E+01 -1.0087E+00 -4.6666E+00  1.1791E+01  6.4585E+00  2.9659E-01  1.1706E+00  8.7866E-01  3.7645E+00  5.1274E+00
            -7.8553E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1603.23861760406        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.9052E-01  6.8233E-01  9.1696E-01  1.2371E+00  7.6214E-01  9.6545E-01  1.5196E+00  2.3101E-01  1.0889E+00  9.0418E-01
             1.1452E+00
 PARAMETER:  9.0478E-02 -2.8225E-01  1.3312E-02  3.1280E-01 -1.7163E-01  6.4843E-02  5.1844E-01 -1.3653E+00  1.8519E-01 -7.2653E-04
             2.3560E-01
 GRADIENT:   7.0442E+00 -1.9768E+00 -2.7937E+00  4.0549E+00  3.7142E+00  3.1677E-02  8.1143E-01  8.2871E-01  1.0299E+00  3.6425E+00
            -6.3435E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1603.23933447132        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.9023E-01  6.8752E-01  9.1091E-01  1.2328E+00  7.6060E-01  9.6554E-01  1.5170E+00  2.2533E-01  1.0894E+00  8.9717E-01
             1.1457E+00
 PARAMETER:  9.0186E-02 -2.7466E-01  6.6916E-03  3.0926E-01 -1.7365E-01  6.4929E-02  5.1672E-01 -1.3902E+00  1.8559E-01 -8.5137E-03
             2.3598E-01
 GRADIENT:   6.0848E+00 -2.2008E+00 -2.6311E+00  2.9691E+00  3.4826E+00 -3.1628E-03  7.3411E-01  7.9629E-01  7.6851E-01  3.3214E+00
            -5.3848E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1603.89010845200        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  1.0013E+00  6.0539E-01  8.8395E-01  1.2906E+00  7.2114E-01  9.6874E-01  1.7213E+00  1.1391E-01  1.0533E+00  8.6137E-01
             1.1514E+00
 PARAMETER:  1.0128E-01 -4.0188E-01 -2.3351E-02  3.5511E-01 -2.2692E-01  6.8238E-02  6.4306E-01 -2.0724E+00  1.5194E-01 -4.9234E-02
             2.4099E-01
 GRADIENT:   3.6382E+00 -1.4335E+00 -3.5009E+00 -1.4996E+00  3.6080E+00 -5.0178E-01 -3.1168E-01  1.8683E-01  1.6674E+00  1.4788E-01
             1.1464E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1604.00769036630        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  9.9868E-01  5.3634E-01  8.9300E-01  1.3323E+00  7.0510E-01  9.6851E-01  1.8835E+00  3.8082E-02  1.0202E+00  8.7467E-01
             1.1478E+00
 PARAMETER:  9.8675E-02 -5.2299E-01 -1.3173E-02  3.8692E-01 -2.4941E-01  6.8009E-02  7.3311E-01 -3.1680E+00  1.2001E-01 -3.3913E-02
             2.3789E-01
 GRADIENT:  -2.2763E-01 -3.5935E-01 -1.7842E-01 -7.0448E-01  2.9102E-01 -3.7827E-03  8.1009E-03  1.8253E-02 -1.3095E-01 -9.4365E-02
             4.4718E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1604.01838476829        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      905
 NPARAMETR:  9.9895E-01  5.5087E-01  8.9482E-01  1.3251E+00  7.1020E-01  9.6895E-01  1.8437E+00  1.2674E-02  1.0264E+00  8.7674E-01
             1.1475E+00
 PARAMETER:  9.8950E-02 -4.9627E-01 -1.1132E-02  3.8148E-01 -2.4220E-01  6.8458E-02  7.1178E-01 -4.2682E+00  1.2605E-01 -3.1540E-02
             2.3758E-01
 GRADIENT:   1.9607E-02 -1.1947E-01 -3.2146E-01  2.2700E-01  6.0010E-01  5.9623E-02 -5.6837E-02  2.0211E-03 -2.0241E-02  2.6899E-02
            -1.8465E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1604.01917369012        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  9.9897E-01  5.5321E-01  8.9263E-01  1.3232E+00  7.0952E-01  9.6888E-01  1.8423E+00  1.0000E-02  1.0268E+00  8.7407E-01
             1.1480E+00
 PARAMETER:  9.9007E-02 -4.9195E-01 -1.3602E-02  3.8011E-01 -2.4316E-01  6.8360E-02  7.1121E-01 -4.6078E+00  1.2653E-01 -3.4580E-02
             2.3805E-01
 GRADIENT:   1.8705E-02  3.5795E-03 -4.9712E-03  2.5739E-02  2.9612E-03 -2.6114E-03  2.9136E-03  0.0000E+00  7.7323E-03  6.0189E-04
             6.6952E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1074
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6392E-05  1.1900E-02 -3.9610E-04 -1.0273E-02 -8.4764E-03
 SE:             2.9773E-02  1.7110E-02  1.9465E-04  2.6625E-02  2.3488E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9849E-01  4.8674E-01  4.1853E-02  6.9963E-01  7.1819E-01

 ETASHRINKSD(%)  2.5794E-01  4.2680E+01  9.9348E+01  1.0802E+01  2.1313E+01
 ETASHRINKVR(%)  5.1522E-01  6.7144E+01  9.9996E+01  2.0438E+01  3.8083E+01
 EBVSHRINKSD(%)  6.0251E-01  4.6061E+01  9.9325E+01  9.8943E+00  1.8476E+01
 EBVSHRINKVR(%)  1.2014E+00  7.0906E+01  9.9995E+01  1.8810E+01  3.3538E+01
 RELATIVEINF(%)  9.8096E+01  3.0617E+00  6.2290E-04  1.3651E+01  5.4348E+00
 EPSSHRINKSD(%)  4.2854E+01
 EPSSHRINKVR(%)  6.7343E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1604.0191736901204     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -868.86834712638222     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1604.019       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  5.53E-01  8.93E-01  1.32E+00  7.10E-01  9.69E-01  1.84E+00  1.00E-02  1.03E+00  8.74E-01  1.15E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.17E+03
 
 TH 2
+       -1.64E+01  3.40E+02
 
 TH 3
+        1.78E+01  2.34E+02  6.58E+02
 
 TH 4
+       -6.53E+00  2.32E+02 -1.27E+02  4.92E+02
 
 TH 5
+       -4.48E+00 -4.90E+02 -1.00E+03  1.37E+02  1.90E+03
 
 TH 6
+        2.61E+00 -4.04E+00  5.72E+00 -1.84E+00 -9.38E-01  2.09E+02
 
 TH 7
+        7.37E-01  2.19E+01  7.48E+00 -1.09E+00 -1.37E+01  1.78E-01  9.02E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.41E-01 -2.27E+01 -1.61E+01  1.17E+01 -1.85E+00  1.83E-01  7.83E+00  0.00E+00  1.29E+02
 
 TH10
+       -3.19E+00  1.96E+01 -7.27E+01 -2.34E+01 -1.77E+01  1.13E+00  9.71E+00  0.00E+00  1.63E+00  1.04E+02
 
 TH11
+       -9.13E+00 -7.93E+00 -3.75E+01 -7.89E+00  1.79E+01  3.55E+00  2.04E+00  0.00E+00  6.06E+00  2.55E+01  1.67E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       17.946
Stop Time:
Sat Sep 18 13:09:46 CDT 2021
