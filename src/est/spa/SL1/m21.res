Sat Sep 18 11:33:29 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat21.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1743.28035128678        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.9722E+01 -3.0367E+01 -1.6415E+01 -3.9386E+01 -4.2232E+01  3.8977E+00  2.2483E+00  1.7605E+01  1.8209E+01  1.5726E+01
             3.1557E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1753.82665097025        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0378E+00  1.0606E+00  1.1768E+00  9.8826E-01  1.1764E+00  9.7671E-01  1.0198E+00  8.7147E-01  8.9907E-01  1.0013E+00
             9.4361E-01
 PARAMETER:  1.3707E-01  1.5884E-01  2.6283E-01  8.8194E-02  2.6246E-01  7.6437E-02  1.1962E-01 -3.7577E-02 -6.3928E-03  1.0127E-01
             4.1956E-02
 GRADIENT:   5.9446E+01 -2.0142E+01 -1.5323E+01 -9.4587E+00  3.5362E+01 -1.2364E+00 -6.2260E+00  4.6411E+00 -1.9396E+00 -1.4946E+01
             6.9583E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1756.12207087165        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0291E+00  9.2913E-01  1.1894E+00  1.0819E+00  1.1087E+00  1.0030E+00  1.2614E+00  6.0849E-01  7.8509E-01  1.0638E+00
             9.3315E-01
 PARAMETER:  1.2868E-01  2.6489E-02  2.7349E-01  1.7875E-01  2.0316E-01  1.0303E-01  3.3220E-01 -3.9677E-01 -1.4195E-01  1.6184E-01
             3.0806E-02
 GRADIENT:   4.0695E+01  2.3165E+00 -2.3356E+00  1.7731E+01  7.2184E+00  9.7254E+00 -1.6128E+00  1.2043E+00 -4.3644E+00 -1.8195E+00
            -2.3113E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1756.74895149297        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0440E+00  8.8549E-01  1.1314E+00  1.1097E+00  1.0530E+00  9.9230E-01  1.3222E+00  4.3031E-01  7.9716E-01  1.0427E+00
             9.3811E-01
 PARAMETER:  1.4305E-01 -2.1616E-02  2.2347E-01  2.0412E-01  1.5162E-01  9.2265E-02  3.7927E-01 -7.4324E-01 -1.2669E-01  1.4177E-01
             3.6114E-02
 GRADIENT:   3.7594E+00  3.8477E+00 -3.4835E-01  4.5321E+00 -2.4171E+00 -1.4324E-01 -7.1283E-02  5.6697E-01  1.6851E-01  6.2365E-01
             4.9092E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1757.01145936747        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      473
 NPARAMETR:  1.0408E+00  7.2743E-01  1.1268E+00  1.1986E+00  9.8795E-01  9.9085E-01  1.5491E+00  1.9156E-01  7.4822E-01  1.0244E+00
             9.3909E-01
 PARAMETER:  1.3999E-01 -2.1824E-01  2.1934E-01  2.8116E-01  8.7881E-02  9.0807E-02  5.3770E-01 -1.5526E+00 -1.9006E-01  1.2414E-01
             3.7158E-02
 GRADIENT:   1.3511E-01 -6.0041E-01  1.0639E-01 -1.8600E+00  3.3162E-02 -6.8877E-02  1.7143E-02  1.7112E-02  1.3437E-01  1.3705E-01
             6.5600E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1757.01416305965        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  1.0405E+00  7.0914E-01  1.1265E+00  1.2101E+00  9.8039E-01  9.9080E-01  1.5804E+00  1.7155E-01  7.4284E-01  1.0217E+00
             9.3913E-01
 PARAMETER:  1.3974E-01 -2.4370E-01  2.1914E-01  2.9066E-01  8.0192E-02  9.0756E-02  5.5768E-01 -1.6629E+00 -1.9727E-01  1.2146E-01
             3.7198E-02
 GRADIENT:  -2.4278E-04 -9.4241E-02 -3.6409E-02 -1.1683E-01  3.0082E-02 -4.2027E-03  4.0681E-03  9.6049E-03  1.3057E-02  1.3665E-02
             1.0774E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1757.01764912755        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  1.0409E+00  7.2761E-01  1.1106E+00  1.1983E+00  9.7970E-01  9.9114E-01  1.5514E+00  7.9052E-02  7.4716E-01  1.0190E+00
             9.3928E-01
 PARAMETER:  1.4004E-01 -2.1799E-01  2.0491E-01  2.8092E-01  7.9494E-02  9.1100E-02  5.3918E-01 -2.4377E+00 -1.9147E-01  1.1880E-01
             3.7360E-02
 GRADIENT:  -2.4031E-02  6.3640E-02 -9.7593E-02  1.6812E-01 -9.8729E-03 -9.7589E-04  2.2013E-03  1.6452E-03  2.1676E-03  4.5432E-02
             3.6750E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1757.01817906711        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1006
 NPARAMETR:  1.0409E+00  7.2823E-01  1.1082E+00  1.1977E+00  9.7882E-01  9.9117E-01  1.5506E+00  1.7474E-02  7.4728E-01  1.0184E+00
             9.3935E-01
 PARAMETER:  1.4008E-01 -2.1713E-01  2.0272E-01  2.8044E-01  7.8593E-02  9.1129E-02  5.3867E-01 -3.9470E+00 -1.9131E-01  1.1825E-01
             3.7429E-02
 GRADIENT:   5.9521E-03  2.3414E-02 -9.9712E-02  1.3861E-01  7.3735E-02  1.2141E-03  7.5330E-05  5.5075E-05  3.3484E-03  1.1668E-02
             1.7875E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1757.01820123617        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1139
 NPARAMETR:  1.0410E+00  7.2752E-01  1.1089E+00  1.1981E+00  9.7879E-01  9.9118E-01  1.5519E+00  1.0000E-02  7.4708E-01  1.0185E+00
             9.3931E-01
 PARAMETER:  1.4006E-01 -2.1819E-01  2.0327E-01  2.8080E-01  7.8578E-02  9.1111E-02  5.3940E-01 -4.7454E+00 -1.9154E-01  1.1847E-01
             3.7450E-02
 GRADIENT:  -3.7668E-02 -5.3096E-03 -7.6047E-03  1.9487E-02  3.9307E-03 -2.1959E-03 -1.5463E-03  0.0000E+00  1.3235E-03  4.0025E-03
             5.0154E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1139
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.4718E-04  9.3762E-03 -3.6457E-04 -1.4080E-02 -1.7350E-02
 SE:             2.9845E-02  1.9991E-02  1.6982E-04  2.3433E-02  2.3800E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9607E-01  6.3905E-01  3.1809E-02  5.4793E-01  4.6600E-01

 ETASHRINKSD(%)  1.6719E-02  3.3028E+01  9.9431E+01  2.1496E+01  2.0269E+01
 ETASHRINKVR(%)  3.3434E-02  5.5147E+01  9.9997E+01  3.8372E+01  3.6429E+01
 EBVSHRINKSD(%)  3.8918E-01  3.4013E+01  9.9459E+01  2.0833E+01  1.7118E+01
 EBVSHRINKVR(%)  7.7685E-01  5.6457E+01  9.9997E+01  3.7326E+01  3.1305E+01
 RELATIVEINF(%)  9.8331E+01  2.6765E+00  3.8662E-04  4.1843E+00  8.4307E+00
 EPSSHRINKSD(%)  4.2245E+01
 EPSSHRINKVR(%)  6.6644E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1757.0182012361668     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1021.8673746724286     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1757.018       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  7.27E-01  1.11E+00  1.20E+00  9.79E-01  9.91E-01  1.55E+00  1.00E-02  7.47E-01  1.02E+00  9.39E-01
 


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
+        1.03E+03
 
 TH 2
+       -1.10E+01  4.00E+02
 
 TH 3
+        1.17E+01  1.04E+02  2.86E+02
 
 TH 4
+       -6.65E+00  4.32E+02 -1.53E+02  8.75E+02
 
 TH 5
+       -2.28E+00 -2.19E+02 -4.03E+02  1.76E+02  7.69E+02
 
 TH 6
+       -1.00E+00 -1.46E+00  3.83E+00 -9.11E-01 -2.09E+00  2.05E+02
 
 TH 7
+        1.05E+00  2.86E+01  5.77E+00 -5.61E+00 -8.33E+00  1.32E-01  2.19E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.55E+00 -1.59E+01 -2.23E+01 -1.27E+01  2.58E+01 -1.25E+00  2.57E+01  0.00E+00  1.48E+02
 
 TH10
+        3.82E-01 -1.83E+00 -2.94E+01 -1.41E+01 -5.55E+01 -6.56E-01  9.60E-01  0.00E+00  2.66E+00  9.62E+01
 
 TH11
+       -7.36E+00 -1.33E+01 -4.15E+01 -6.13E+00  1.44E+01  5.39E+00  2.56E+00  0.00E+00  8.80E+00  2.21E+01  2.64E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.630
Stop Time:
Sat Sep 18 11:33:51 CDT 2021
