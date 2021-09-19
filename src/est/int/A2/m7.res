Sat Sep 18 00:27:53 CDT 2021
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
$DATA ../../../../data/int/A2/dat7.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2192.09750465768        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0099E+02  1.0977E+02  1.7783E+02  7.3749E+01  1.7297E+02 -3.2989E+01 -1.6325E+02 -3.3567E+02 -9.9771E+00 -1.5036E+02
            -2.6100E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3026.56040022403        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8199E-01  1.0670E+00  9.5009E-01  9.5831E-01  9.9753E-01  1.1752E+00  9.4608E-01  1.0025E+00  7.4187E-01  1.0863E+00
             1.9842E+00
 PARAMETER:  8.1830E-02  1.6484E-01  4.8803E-02  5.7419E-02  9.7530E-02  2.6140E-01  4.4572E-02  1.0249E-01 -1.9858E-01  1.8278E-01
             7.8522E-01
 GRADIENT:   1.0672E+01  1.9873E+01 -1.0257E+01  1.5801E+01  3.2641E+00  3.4741E+01 -6.7059E+00  2.7900E+00 -2.3584E+01 -1.3526E+01
            -7.7905E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3032.34066027186        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.8387E-01  1.0801E+00  8.9554E-01  9.4630E-01  9.9077E-01  1.1162E+00  1.0164E+00  2.8726E-01  7.8309E-01  1.1527E+00
             2.0584E+00
 PARAMETER:  8.3740E-02  1.7704E-01 -1.0331E-02  4.4806E-02  9.0731E-02  2.0997E-01  1.1624E-01 -1.1474E+00 -1.4451E-01  2.4208E-01
             8.2194E-01
 GRADIENT:   1.2486E+01  2.2963E+01 -1.2145E+01  2.6148E+00  6.6341E+00  1.7400E+01  6.9429E+00  6.1697E-01 -1.1393E+01 -5.1245E+00
            -4.4961E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3034.75909542123        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.7842E-01  9.7570E-01  8.7699E-01  9.9558E-01  9.1080E-01  1.0663E+00  9.0867E-01  1.6263E-01  8.5966E-01  1.1395E+00
             2.0608E+00
 PARAMETER:  7.8181E-02  7.5397E-02 -3.1255E-02  9.5573E-02  6.5730E-03  1.6421E-01  4.2247E-03 -1.7162E+00 -5.1220E-02  2.3060E-01
             8.2311E-01
 GRADIENT:   2.1298E+00 -2.8150E+00  1.9806E-02 -8.9039E-01  2.4652E-01  2.2248E-01  2.2399E-01  2.9987E-01  1.2469E+00 -2.6210E-02
             1.6597E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3034.91965878926        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.7749E-01  9.8592E-01  8.8053E-01  9.9119E-01  9.1901E-01  1.0657E+00  9.0293E-01  3.9522E-02  8.5835E-01  1.1492E+00
             2.0595E+00
 PARAMETER:  7.7238E-02  8.5815E-02 -2.7233E-02  9.1147E-02  1.5541E-02  1.6365E-01 -2.1157E-03 -3.1309E+00 -5.2740E-02  2.3904E-01
             8.2248E-01
 GRADIENT:   3.1333E-01 -1.0568E-01 -2.1957E-01  3.0706E-01  6.2669E-01 -3.9558E-03  1.3182E-01  1.7448E-02  4.2521E-01  2.6169E-01
            -5.4395E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3034.92024747550        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.7738E-01  9.8443E-01  8.7969E-01  9.9191E-01  9.1736E-01  1.0657E+00  9.0350E-01  2.5486E-02  8.5695E-01  1.1473E+00
             2.0599E+00
 PARAMETER:  7.7123E-02  8.4305E-02 -2.8182E-02  9.1874E-02  1.3742E-02  1.6364E-01 -1.4814E-03 -3.5696E+00 -5.4377E-02  2.3744E-01
             8.2265E-01
 GRADIENT:   8.0565E-02 -2.2406E-02 -3.2206E-02  9.7166E-02  1.4061E-01 -8.5196E-03  4.4487E-02  7.3258E-03  1.1290E-01  6.7714E-02
            -1.6935E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3034.92109861895        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  9.7734E-01  9.8393E-01  8.7934E-01  9.9213E-01  9.1679E-01  1.0657E+00  9.0355E-01  1.2858E-02  8.5649E-01  1.1467E+00
             2.0600E+00
 PARAMETER:  7.7084E-02  8.3801E-02 -2.8585E-02  9.2097E-02  1.3126E-02  1.6366E-01 -1.4231E-03 -4.2538E+00 -5.4914E-02  2.3692E-01
             8.2273E-01
 GRADIENT:   1.7826E-03 -9.4239E-03  2.7066E-03 -8.4247E-03 -5.2238E-03 -1.4835E-03  6.6751E-03  1.8855E-03  2.7044E-03 -5.0868E-04
            -7.4019E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3034.96927207537        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      577
 NPARAMETR:  9.8299E-01  9.9509E-01  8.8707E-01  9.8770E-01  9.2765E-01  1.0710E+00  9.0124E-01  1.0000E-02  8.5705E-01  1.1543E+00
             2.0622E+00
 PARAMETER:  8.2844E-02  9.5082E-02 -1.9834E-02  8.7624E-02  2.4897E-02  1.6860E-01 -3.9844E-03 -4.9007E+00 -5.4262E-02  2.4353E-01
             8.2378E-01
 GRADIENT:   7.2701E-01  2.0036E-02  2.0843E-01 -2.7160E-01 -2.7073E-01  1.5437E-01  2.6945E-02  0.0000E+00 -1.8881E-01 -2.1369E-01
             3.6424E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -3034.96983398269        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  9.8261E-01  9.9577E-01  8.8716E-01  9.8741E-01  9.2833E-01  1.0706E+00  8.9931E-01  1.0000E-02  8.5815E-01  1.1562E+00
             2.0619E+00
 PARAMETER:  8.2460E-02  9.5762E-02 -1.9730E-02  8.7332E-02  2.5628E-02  1.6820E-01 -6.1257E-03 -4.9044E+00 -5.2974E-02  2.4513E-01
             8.2361E-01
 GRADIENT:  -1.2354E-03 -1.8445E-03  3.0307E-03 -2.9437E-03 -3.2583E-04 -1.2313E-03  1.5106E-04  0.0000E+00  8.4575E-04 -1.2173E-03
            -3.4804E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      669
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.9282E-04 -1.7255E-02 -9.0354E-05  5.5052E-03 -1.1799E-02
 SE:             2.9635E-02  1.9591E-02  1.3950E-04  2.5674E-02  2.6121E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7597E-01  3.7844E-01  5.1716E-01  8.3021E-01  6.5147E-01

 ETASHRINKSD(%)  7.1973E-01  3.4367E+01  9.9533E+01  1.3990E+01  1.2492E+01
 ETASHRINKVR(%)  1.4343E+00  5.6923E+01  9.9998E+01  2.6023E+01  2.3424E+01
 EBVSHRINKSD(%)  8.9741E-01  3.4430E+01  9.9496E+01  1.4561E+01  1.2784E+01
 EBVSHRINKVR(%)  1.7868E+00  5.7005E+01  9.9997E+01  2.7001E+01  2.3934E+01
 RELATIVEINF(%)  9.8189E+01  9.5352E+00  1.3092E-03  2.4592E+01  1.5738E+01
 EPSSHRINKSD(%)  1.7217E+01
 EPSSHRINKVR(%)  3.1470E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3034.9698339826891     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1380.8804742142784     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.01
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3034.970       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  9.96E-01  8.87E-01  9.87E-01  9.28E-01  1.07E+00  8.99E-01  1.00E-02  8.58E-01  1.16E+00  2.06E+00
 


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
+        9.85E+02
 
 TH 2
+       -6.23E+00  5.97E+02
 
 TH 3
+       -1.02E-02  5.20E+01  3.29E+02
 
 TH 4
+       -1.48E+01  4.36E+02 -9.82E+01  1.12E+03
 
 TH 5
+       -9.80E-01 -3.72E+02 -2.97E+02  1.84E+02  7.26E+02
 
 TH 6
+        2.14E+00 -1.55E+00  2.40E+00 -5.41E+00 -9.81E-01  1.65E+02
 
 TH 7
+       -7.69E-01  3.77E+00 -5.23E+00  2.65E+00 -3.10E+00 -4.48E-02  3.98E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.09E+00 -1.44E+01  7.62E+00  2.76E+01  2.61E+00 -2.66E-01  2.32E+01  0.00E+00  1.43E+02
 
 TH10
+        3.75E-01 -1.61E+01 -1.71E+01  5.08E+00 -1.73E+01  2.85E-01  2.18E+01  0.00E+00  7.50E+00  8.21E+01
 
 TH11
+       -1.15E+01 -1.72E+01 -5.40E+00 -1.76E+01  2.50E+00  1.81E+00  5.22E+00  0.00E+00  6.90E+00  5.01E+00  2.76E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.250
Stop Time:
Sat Sep 18 00:28:19 CDT 2021
