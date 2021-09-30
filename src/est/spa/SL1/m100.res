Wed Sep 29 15:31:25 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat100.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1664.39373382677        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0559E+02 -1.8978E+01 -2.8218E+01  5.2242E+01  4.5457E+01  3.2797E+01  1.1495E+00 -1.0262E+00  3.7850E+01  7.3909E+00
             1.2206E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1671.73739532106        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.6283E-01  1.0392E+00  9.5874E-01  9.9810E-01  9.5832E-01  1.0236E+00  1.0932E+00  1.0570E+00  7.9661E-01  9.3206E-01
             1.0241E+00
 PARAMETER:  6.2123E-02  1.3840E-01  5.7863E-02  9.8097E-02  5.7422E-02  1.2332E-01  1.8911E-01  1.5539E-01 -1.2739E-01  2.9645E-02
             1.2378E-01
 GRADIENT:   3.9799E+02  4.0727E+01 -6.1456E+00  8.2427E+01 -8.9102E-01  6.0631E+01 -2.1639E-01  2.1105E+00  1.0314E+01  7.0450E+00
             2.1461E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1673.16419415367        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      248
 NPARAMETR:  9.6368E-01  1.0848E+00  9.8580E-01  9.6325E-01  1.0050E+00  1.0458E+00  1.1975E+00  1.0710E+00  7.8566E-01  9.4408E-01
             9.6234E-01
 PARAMETER:  6.3002E-02  1.8141E-01  8.5698E-02  6.2558E-02  1.0501E-01  1.4475E-01  2.8022E-01  1.6863E-01 -1.4123E-01  4.2458E-02
             6.1610E-02
 GRADIENT:   2.1477E+00 -8.2924E+00 -5.0625E+00 -5.6530E-01  6.0946E+00  3.8653E+00  6.2193E+00 -1.1328E+00  8.5437E+00  4.1916E-01
            -5.1353E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1674.30851573544        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  9.6353E-01  1.1719E+00  1.2377E+00  9.2250E-01  1.1349E+00  1.0330E+00  1.1474E+00  1.5549E+00  7.1131E-01  1.0464E+00
             9.7712E-01
 PARAMETER:  6.2848E-02  2.5863E-01  3.1323E-01  1.9335E-02  2.2659E-01  1.3246E-01  2.3746E-01  5.4138E-01 -2.4065E-01  1.4540E-01
             7.6852E-02
 GRADIENT:   3.0010E+00  3.0088E+00  1.3819E+00  2.7315E+00 -4.1779E+00 -9.0633E-01  1.9601E+00  3.9007E-01 -1.6291E-01  7.6193E-02
            -7.4935E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1674.58120597739        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      601
 NPARAMETR:  9.6206E-01  1.3723E+00  1.0280E+00  7.9535E-01  1.1530E+00  1.0426E+00  9.7329E-01  1.4989E+00  8.1520E-01  1.0455E+00
             9.8000E-01
 PARAMETER:  6.1326E-02  4.1649E-01  1.2757E-01 -1.2897E-01  2.4235E-01  1.4174E-01  7.2930E-02  5.0471E-01 -1.0432E-01  1.4447E-01
             7.9798E-02
 GRADIENT:  -3.0752E+00  8.2012E+00  2.3154E+00  7.2188E+00 -3.8785E+00  2.1276E+00 -2.2398E+00 -8.6717E-01 -4.9173E-01  2.3673E-01
             2.0512E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1674.73625851689        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      778
 NPARAMETR:  9.6500E-01  1.5200E+00  8.0404E-01  6.9269E-01  1.1409E+00  1.0346E+00  9.1047E-01  1.3875E+00  8.8017E-01  1.0166E+00
             9.7920E-01
 PARAMETER:  6.4372E-02  5.1869E-01 -1.1811E-01 -2.6717E-01  2.3185E-01  1.3402E-01  6.2079E-03  4.2748E-01 -2.7641E-02  1.1649E-01
             7.8979E-02
 GRADIENT:   7.6227E-01  4.4496E+00  1.3709E+00  2.7738E+00 -3.0460E+00 -1.4001E+00  1.0478E-01 -5.5008E-01  5.4262E-02 -3.1837E-01
             6.0010E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1674.73857999284        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      954
 NPARAMETR:  9.6546E-01  1.5795E+00  7.3618E-01  6.5362E-01  1.1439E+00  1.0342E+00  8.8427E-01  1.3693E+00  9.0921E-01  1.0134E+00
             9.7927E-01
 PARAMETER:  6.4848E-02  5.5712E-01 -2.0628E-01 -3.2523E-01  2.3446E-01  1.3364E-01 -2.2995E-02  4.1428E-01  4.8197E-03  1.1329E-01
             7.9051E-02
 GRADIENT:   1.0162E+00  5.6948E+00  1.1815E+00  3.4651E+00 -2.8569E+00 -1.7069E+00  1.0489E-01 -5.3267E-01 -1.2623E-01 -4.2645E-01
             1.0872E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1674.73931158351        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1139
 NPARAMETR:  9.6554E-01  1.6057E+00  7.1034E-01  6.3632E-01  1.1475E+00  1.0345E+00  8.7295E-01  1.3697E+00  9.2365E-01  1.0143E+00
             9.7934E-01
 PARAMETER:  6.4937E-02  5.7356E-01 -2.4201E-01 -3.5206E-01  2.3758E-01  1.3395E-01 -3.5878E-02  4.1459E-01  2.0577E-02  1.1417E-01
             7.9123E-02
 GRADIENT:   9.5152E-01  5.6838E+00  1.0900E+00  3.4962E+00 -2.6548E+00 -1.6381E+00  6.9080E-02 -5.1733E-01 -1.5844E-01 -4.2336E-01
            -2.3516E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1674.75493637913        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1303
 NPARAMETR:  9.6604E-01  1.6120E+00  7.0794E-01  6.2921E-01  1.1491E+00  1.0361E+00  8.7044E-01  1.4191E+00  9.2504E-01  1.0149E+00
             9.7885E-01
 PARAMETER:  6.5445E-02  5.7745E-01 -2.4540E-01 -3.6330E-01  2.3899E-01  1.3543E-01 -3.8752E-02  4.5005E-01  2.2080E-02  1.1483E-01
             7.8625E-02
 GRADIENT:   2.0730E+00  2.4333E-01  1.4423E+00 -5.3966E-01 -5.4251E+00 -1.0250E+00  3.3839E-01  2.3310E-02 -8.1858E-02 -2.7378E-02
            -6.9303E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1674.77149580724        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1438
 NPARAMETR:  9.6561E-01  1.6096E+00  7.0794E-01  6.2978E-01  1.1532E+00  1.0392E+00  8.6935E-01  1.4192E+00  9.2515E-01  1.0186E+00
             9.7896E-01
 PARAMETER:  6.5002E-02  5.7600E-01 -2.4540E-01 -3.6238E-01  2.4257E-01  1.3845E-01 -4.0012E-02  4.5006E-01  2.2197E-02  1.1840E-01
             7.8733E-02
 GRADIENT:  -6.0618E-01  6.0399E-01  3.2231E+04  2.1791E+04  3.2581E+04 -2.2473E-01 -5.6708E-02 -2.6775E+02  2.5274E-02  6.6745E+04
             9.6710E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1438
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.3149E-04 -1.4539E-02 -3.4385E-02  1.2182E-02 -3.7832E-02
 SE:             2.9871E-02  2.5088E-02  1.2563E-02  1.8699E-02  2.1412E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8580E-01  5.6225E-01  6.1984E-03  5.1474E-01  7.7251E-02

 ETASHRINKSD(%)  1.0000E-10  1.5950E+01  5.7913E+01  3.7357E+01  2.8267E+01
 ETASHRINKVR(%)  1.0000E-10  2.9357E+01  8.2287E+01  6.0758E+01  4.8544E+01
 EBVSHRINKSD(%)  3.8571E-01  1.5856E+01  6.1543E+01  4.0092E+01  2.4817E+01
 EBVSHRINKVR(%)  7.6993E-01  2.9198E+01  8.5211E+01  6.4110E+01  4.3475E+01
 RELATIVEINF(%)  9.8964E+01  2.0054E+00  6.7946E-01  8.1893E-01  1.3105E+01
 EPSSHRINKSD(%)  4.4830E+01
 EPSSHRINKVR(%)  6.9563E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1674.7714958072370     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.62066924349881     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1674.771       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.66E-01  1.61E+00  7.08E-01  6.30E-01  1.15E+00  1.04E+00  8.69E-01  1.42E+00  9.25E-01  1.02E+00  9.79E-01
 


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
+        1.10E+03
 
 TH 2
+       -6.74E+00  2.30E+05
 
 TH 3
+        1.18E+07  6.21E+01  6.53E+06
 
 TH 4
+        1.11E+02  9.34E+05 -8.87E+03  7.59E+06
 
 TH 5
+        7.31E+06 -7.62E+05 -7.25E+03  6.15E+02  2.53E+06
 
 TH 6
+        6.88E-01 -1.74E+00  1.61E+00  1.18E+02  9.84E+01  1.82E+02
 
 TH 7
+        1.66E+00  1.72E+01 -1.31E+07  7.67E+02  6.40E+02  2.08E-01  1.38E+02
 
 TH 8
+       -3.18E+06 -8.51E+00 -1.77E+06  2.07E+04 -1.10E+06 -2.13E+06  2.89E+00  4.77E+05
 
 TH 9
+        6.66E-01 -1.05E+01  1.23E+07 -2.68E+03  7.63E+06 -7.38E-02  3.40E+01  3.97E+01  3.47E+01
 
 TH10
+        1.70E+07 -1.77E+06 -1.63E+04  7.17E+06  5.86E+06  2.31E+02  1.49E+03  3.91E+04 -5.12E+03  1.36E+07
 
 TH11
+       -2.09E+07  2.18E+06 -1.16E+07 -3.03E+03 -7.22E+06  1.56E+00  8.01E+00  3.14E+06  8.06E+00 -1.67E+07  2.06E+07
 
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
 #CPUT: Total CPU Time in Seconds,       25.081
Stop Time:
Wed Sep 29 15:31:52 CDT 2021
