Sat Sep 18 02:54:52 CDT 2021
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
$DATA ../../../../data/int/S1/dat81.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3217.63773584551        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.9787E+00 -1.4480E+01  8.3887E+01 -2.0093E+01  1.0411E+02 -3.8592E+01 -4.2910E+01 -2.5720E+02 -3.6215E+01 -4.8751E+01
            -9.7255E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3759.09729404868        NO. OF FUNC. EVALS.:  87
 CUMULATIVE NO. OF FUNC. EVALS.:      100
 NPARAMETR:  9.8585E-01  1.0167E+00  9.3402E-01  9.9919E-01  9.6588E-01  1.1366E+00  9.8781E-01  1.1340E+00  9.9347E-01  1.0646E+00
             1.2813E+00
 PARAMETER:  8.5747E-02  1.1658E-01  3.1743E-02  9.9187E-02  6.5280E-02  2.2808E-01  8.7737E-02  2.2576E-01  9.3453E-02  1.6265E-01
             3.4787E-01
 GRADIENT:  -2.9911E+01 -2.5315E+01 -2.8372E+01  3.0380E+00  2.1832E+01  1.9051E+01  7.6077E+00  6.3438E+00  8.7927E+00  3.8092E+00
             2.7581E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3759.76045703055        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      281
 NPARAMETR:  9.8587E-01  1.0167E+00  9.3404E-01  9.9921E-01  9.6590E-01  1.1212E+00  9.3093E-01  1.0189E+00  9.6487E-01  1.0646E+00
             1.2814E+00
 PARAMETER:  8.5772E-02  1.1655E-01  3.1768E-02  9.9212E-02  6.5305E-02  2.1436E-01  2.8433E-02  1.1872E-01  6.4241E-02  1.6261E-01
             3.4795E-01
 GRADIENT:  -6.1554E+01 -3.4986E+01 -2.7157E+01 -5.4200E+00  2.2799E+01 -2.5588E-01 -8.5542E-02 -2.4209E-02  2.3930E-02  1.9750E-01
             2.6559E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3759.76433590014        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  9.8587E-01  1.0167E+00  9.3404E-01  9.9921E-01  9.6590E-01  1.1228E+00  9.2592E-01  1.0139E+00  9.6676E-01  1.0646E+00
             1.2814E+00
 PARAMETER:  8.5774E-02  1.1655E-01  3.1769E-02  9.9213E-02  6.5306E-02  2.1578E-01  2.3033E-02  1.1385E-01  6.6199E-02  1.6260E-01
             3.4793E-01
 GRADIENT:  -6.1369E+01 -3.5441E+01 -2.6969E+01 -5.3806E+00  2.3006E+01  3.1210E-01 -6.5432E-01 -2.9225E-01  4.6196E-01 -7.2249E-02
             2.6517E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3771.19250193547        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      645
 NPARAMETR:  9.9218E-01  1.0150E+00  9.3808E-01  1.0022E+00  9.6783E-01  1.1465E+00  6.7598E-01  1.0481E+00  1.0473E+00  1.0599E+00
             1.1743E+00
 PARAMETER:  9.2145E-02  1.1486E-01  3.6080E-02  1.0222E-01  6.7298E-02  2.3671E-01 -2.9159E-01  1.4699E-01  1.4619E-01  1.5821E-01
             2.6067E-01
 GRADIENT:  -4.6664E+01 -5.9482E+01 -1.0507E+01 -4.7423E-01  3.1758E+01  8.9842E+00 -2.6977E+01 -7.6400E+00  1.4427E+01 -2.3577E+01
             1.0933E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3774.91940723352        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      829
 NPARAMETR:  9.9342E-01  1.0146E+00  9.3887E-01  1.0028E+00  9.6821E-01  1.1424E+00  7.1688E-01  1.1030E+00  1.0364E+00  1.0590E+00
             1.1544E+00
 PARAMETER:  9.3395E-02  1.1453E-01  3.6925E-02  1.0280E-01  6.7689E-02  2.3314E-01 -2.3285E-01  1.9807E-01  1.3571E-01  1.5735E-01
             2.4356E-01
 GRADIENT:  -4.4686E+01 -5.2383E+01 -1.1085E+01  3.5210E-01  3.0432E+01  7.6477E+00 -2.4744E+01 -6.0288E+00  1.1963E+01 -2.0287E+01
             8.5713E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3776.36301972176        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      986
 NPARAMETR:  9.9341E-01  1.0320E+00  9.3889E-01  1.0028E+00  9.6819E-01  1.1241E+00  7.1715E-01  1.1031E+00  1.0363E+00  1.1098E+00
             1.1544E+00
 PARAMETER:  9.3387E-02  1.3154E-01  3.6940E-02  1.0279E-01  6.7677E-02  2.1695E-01 -2.3247E-01  1.9815E-01  1.3564E-01  2.0414E-01
             2.4354E-01
 GRADIENT:  -4.6270E+01 -2.5596E+01 -1.0200E+01  1.2996E+01  1.2560E+01  1.3735E+00 -2.0306E+01 -6.3559E+00  1.0733E+01 -9.5649E+00
             8.5221E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3776.37881452328        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1172
 NPARAMETR:  9.9346E-01  1.0320E+00  9.3893E-01  1.0028E+00  9.6815E-01  1.1201E+00  7.1723E-01  1.1030E+00  1.0364E+00  1.1099E+00
             1.1542E+00
 PARAMETER:  9.3437E-02  1.3147E-01  3.6990E-02  1.0284E-01  6.7627E-02  2.1345E-01 -2.3236E-01  1.9805E-01  1.3571E-01  2.0424E-01
             2.4342E-01
 GRADIENT:  -4.6511E+01 -2.5610E+01 -1.0146E+01  1.3046E+01  1.2532E+01 -1.4861E-03 -2.0302E+01 -6.3717E+00  1.0758E+01 -9.5292E+00
             8.4987E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3776.51402134815        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1370
 NPARAMETR:  9.9377E-01  1.0322E+00  9.3914E-01  1.0028E+00  9.6804E-01  1.1201E+00  7.1939E-01  1.1044E+00  1.0359E+00  1.1104E+00
             1.1538E+00
 PARAMETER:  9.3746E-02  1.3166E-01  3.7205E-02  1.0275E-01  6.7515E-02  2.1342E-01 -2.2935E-01  1.9929E-01  1.3525E-01  2.0472E-01
             2.4302E-01
 GRADIENT:  -4.5964E+01 -2.5054E+01 -1.0046E+01  1.2976E+01  1.2043E+01  1.4030E-02 -2.0119E+01 -6.3134E+00  1.0664E+01 -9.2408E+00
             8.4493E+01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -3776.51402134815        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1400
 NPARAMETR:  9.9377E-01  1.0322E+00  9.3914E-01  1.0028E+00  9.6804E-01  1.1201E+00  7.1939E-01  1.1044E+00  1.0359E+00  1.1104E+00
             1.1538E+00
 PARAMETER:  9.3746E-02  1.3166E-01  3.7205E-02  1.0275E-01  6.7515E-02  2.1342E-01 -2.2935E-01  1.9929E-01  1.3525E-01  2.0472E-01
             2.4302E-01
 GRADIENT:  -1.1003E+06 -8.3567E+05 -1.1002E+06  1.0708E+06 -2.2004E+06 -2.7614E-03  9.5940E+05  1.1041E+06  8.1348E+05  1.0748E+06
            -9.0567E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1400
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3594E-02 -2.5222E-02 -1.9637E-02  6.3469E-03 -3.0448E-02
 SE:             2.9809E-02  2.0874E-02  1.9117E-02  2.7003E-02  2.6753E-02
 N:                     100         100         100         100         100

 P VAL.:         4.2866E-01  2.2692E-01  3.0432E-01  8.1417E-01  2.5507E-01

 ETASHRINKSD(%)  1.3468E-01  3.0071E+01  3.5956E+01  9.5370E+00  1.0375E+01
 ETASHRINKVR(%)  2.6917E-01  5.1099E+01  5.8984E+01  1.8165E+01  1.9674E+01
 EBVSHRINKSD(%)  2.7673E-01  3.7978E+01  4.0624E+01  6.6627E+00  1.2735E+01
 EBVSHRINKVR(%)  5.5270E-01  6.1533E+01  6.4745E+01  1.2881E+01  2.3849E+01
 RELATIVEINF(%)  9.9443E+01  1.5088E+01  2.4470E+01  5.4085E+01  2.9546E+01
 EPSSHRINKSD(%)  2.3857E+01
 EPSSHRINKVR(%)  4.2022E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3776.5140213481536     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2122.4246615797429     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.71
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3776.514       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  1.03E+00  9.39E-01  1.00E+00  9.68E-01  1.12E+00  7.19E-01  1.10E+00  1.04E+00  1.11E+00  1.15E+00
 


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
+        5.57E+09
 
 TH 2
+       -3.59E+03  2.98E+09
 
 TH 3
+        5.89E+09 -1.09E+05  6.24E+09
 
 TH 4
+       -5.37E+09  4.51E+04 -5.68E+09  5.18E+09
 
 TH 5
+        5.72E+09 -5.05E+04  6.05E+09 -5.52E+09  5.87E+09
 
 TH 6
+        2.32E+09 -6.00E+03  2.45E+09 -2.23E+09 -8.42E+03  1.57E+02
 
 TH 7
+        2.96E+03  6.78E+04  8.97E+04 -3.68E+04  4.13E+04  4.94E+03  2.02E+09
 
 TH 8
+        2.22E+03  5.08E+04  6.73E+04 -2.76E+04  3.09E+04  3.70E+03 -2.44E+05  1.14E+09
 
 TH 9
+        3.49E+03 -2.89E+09  1.06E+05 -4.34E+04 -4.06E+09  5.82E+03 -2.41E+04 -1.80E+04  2.80E+09
 
 TH10
+        2.15E+03 -1.78E+09  6.52E+04 -2.67E+04 -2.50E+09  3.59E+03 -2.81E+04 -2.11E+04  1.73E+09  1.06E+09
 
 TH11
+       -1.75E+03 -4.00E+04 -5.29E+04  2.17E+04 -2.43E+04 -2.91E+03  7.14E+04  1.43E+05  1.42E+04  1.66E+04  7.00E+08
 
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
 #CPUT: Total CPU Time in Seconds,       52.837
Stop Time:
Sat Sep 18 02:55:46 CDT 2021
