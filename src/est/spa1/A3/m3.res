Wed Sep 29 23:54:18 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   674.654509413139        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8277E+02  4.0607E+01  2.0152E+02 -5.7825E+01  2.7521E+02  5.8847E+01 -5.3116E+01 -3.2851E+02 -1.2011E+02 -1.6454E+02
            -4.8770E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1449.78800183044        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1191E+00  1.0972E+00  9.6819E-01  1.1697E+00  9.3841E-01  8.3718E-01  9.2844E-01  1.0063E+00  8.5930E-01  9.6516E-01
             5.3695E+00
 PARAMETER:  2.1249E-01  1.9272E-01  6.7672E-02  2.5679E-01  3.6428E-02 -7.7715E-02  2.5755E-02  1.0632E-01 -5.1637E-02  6.4544E-02
             1.7807E+00
 GRADIENT:   1.4230E+02  2.3100E+01 -9.3707E+00  4.9659E+01 -2.0512E+01 -9.2538E+00  1.3263E+01  7.1519E+00  2.4263E+01  2.1290E+01
             2.8021E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1492.16901578503        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0499E+00  9.4139E-01  4.6814E-01  1.1211E+00  6.1765E-01  8.5822E-01  3.6301E-01  1.0664E-01  8.4863E-01  6.4042E-01
             4.6472E+00
 PARAMETER:  1.4870E-01  3.9597E-02 -6.5899E-01  2.1427E-01 -3.8183E-01 -5.2892E-02 -9.1334E-01 -2.1383E+00 -6.4126E-02 -3.4563E-01
             1.6363E+00
 GRADIENT:   6.5176E+00 -2.7929E+01 -3.6446E+01  2.3775E+01  4.4857E+01 -9.1116E+00  4.1038E-01  1.8110E-01  1.0383E+01  1.8693E+01
             1.9434E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1518.53156182769        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0270E+00  9.8884E-01  4.9142E-01  1.0564E+00  6.3539E-01  9.0728E-01  6.5153E-01  5.5945E-02  8.6341E-01  2.5862E-01
             3.7701E+00
 PARAMETER:  1.2660E-01  8.8774E-02 -6.1045E-01  1.5487E-01 -3.5351E-01  2.6960E-03 -3.2843E-01 -2.7834E+00 -4.6869E-02 -1.2524E+00
             1.4271E+00
 GRADIENT:   2.4887E+00 -8.4855E+00 -3.1152E+00 -7.9876E+00  3.9544E+00  5.7374E+00 -1.5426E+00  2.4059E-02  2.7551E+00  1.8653E+00
            -8.8881E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1524.14786365452        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      344
 NPARAMETR:  1.0345E+00  9.7806E-01  1.2238E+00  1.1558E+00  1.0700E+00  8.6998E-01  7.4187E-01  7.5466E-02  7.1602E-01  1.7862E-01
             3.9440E+00
 PARAMETER:  1.3390E-01  7.7811E-02  3.0197E-01  2.4482E-01  1.6769E-01 -3.9285E-02 -1.9859E-01 -2.4841E+00 -2.3405E-01 -1.6225E+00
             1.4722E+00
 GRADIENT:   1.8545E-01 -1.1661E+00 -3.0956E+00 -8.6729E-01  5.2706E+00 -1.9077E-01  6.5900E-01  3.1624E-02 -9.4548E-01  5.8165E-01
             4.8534E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1524.50742018762        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      520
 NPARAMETR:  1.0330E+00  8.0958E-01  1.2659E+00  1.2597E+00  1.0094E+00  8.7189E-01  6.5454E-01  7.3095E-02  7.1372E-01  8.1965E-02
             3.9213E+00
 PARAMETER:  1.3243E-01 -1.1124E-01  3.3581E-01  3.3085E-01  1.0935E-01 -3.7094E-02 -3.2382E-01 -2.5160E+00 -2.3727E-01 -2.4015E+00
             1.4664E+00
 GRADIENT:   7.7321E-01  1.9721E+00  4.4653E-01  3.4946E+00 -1.3254E+00  1.8042E-01 -6.8043E-01  3.4564E-02 -3.8872E-01  1.1281E-01
            -1.8143E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1524.68424203478        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      695
 NPARAMETR:  1.0308E+00  6.5405E-01  1.2400E+00  1.3522E+00  9.4130E-01  8.7130E-01  8.9692E-01  4.0768E-02  6.6258E-01  2.5054E-02
             3.9167E+00
 PARAMETER:  1.3029E-01 -3.2457E-01  3.1514E-01  4.0172E-01  3.9503E-02 -3.7764E-02 -8.7922E-03 -3.0999E+00 -3.1162E-01 -3.5867E+00
             1.4652E+00
 GRADIENT:  -3.8043E-01  1.4152E+00 -1.7852E-01  3.9061E+00 -2.1708E-01  1.8816E-01  8.0941E-02  1.2942E-02  2.9047E-02  1.1099E-02
            -1.4108E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1524.70506200831        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      872
 NPARAMETR:  1.0294E+00  5.5553E-01  1.2588E+00  1.4097E+00  9.1628E-01  8.6999E-01  9.8101E-01  2.7510E-02  6.4234E-01  1.0000E-02
             3.9195E+00
 PARAMETER:  1.2900E-01 -4.8783E-01  3.3020E-01  4.4339E-01  1.2564E-02 -3.9270E-02  8.0831E-02 -3.4932E+00 -3.4264E-01 -4.6853E+00
             1.4660E+00
 GRADIENT:   2.9045E-02  1.3690E-01 -1.1581E-02  4.1265E-01 -1.2453E-02  3.7341E-03  1.5372E-03  6.1043E-03 -5.0217E-02  0.0000E+00
            -9.1057E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1524.70777773439        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1005
 NPARAMETR:  1.0294E+00  5.5541E-01  1.2590E+00  1.4096E+00  9.1618E-01  8.6996E-01  9.7985E-01  1.0000E-02  6.4258E-01  1.0000E-02
             3.9195E+00
 PARAMETER:  1.2899E-01 -4.8806E-01  3.3031E-01  4.4329E-01  1.2453E-02 -3.9311E-02  7.9639E-02 -4.5339E+00 -3.4226E-01 -4.6853E+00
             1.4660E+00
 GRADIENT:   4.5515E-02  4.5350E-02  5.6224E-02 -5.3472E-03 -6.4178E-02 -2.0013E-04  1.1940E-03  1.5293E-04 -1.7603E-02  0.0000E+00
            -4.6868E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1005
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0671E-03 -7.2161E-03  6.5358E-05 -1.0619E-02 -1.6314E-05
 SE:             2.8473E-02  8.9707E-03  8.1565E-05  2.2160E-02  1.4132E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4213E-01  4.2116E-01  4.2296E-01  6.3179E-01  9.0810E-01

 ETASHRINKSD(%)  4.6135E+00  6.9947E+01  9.9727E+01  2.5762E+01  9.9527E+01
 ETASHRINKVR(%)  9.0141E+00  9.0968E+01  9.9999E+01  4.4887E+01  9.9998E+01
 EBVSHRINKSD(%)  4.5959E+00  7.0004E+01  9.9661E+01  2.5872E+01  9.9443E+01
 EBVSHRINKVR(%)  8.9807E+00  9.1002E+01  9.9999E+01  4.5050E+01  9.9997E+01
 RELATIVEINF(%)  7.6003E+01  4.5748E-02  4.1710E-05  3.3838E-01  9.2830E-05
 EPSSHRINKSD(%)  1.5308E+01
 EPSSHRINKVR(%)  2.8272E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1524.7077777343925     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -605.76924452971980     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1524.708       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  5.55E-01  1.26E+00  1.41E+00  9.16E-01  8.70E-01  9.80E-01  1.00E-02  6.43E-01  1.00E-02  3.92E+00
 


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
+        1.26E+03
 
 TH 2
+       -8.87E+01  3.35E+02
 
 TH 3
+       -4.21E+00  4.71E+01  7.43E+01
 
 TH 4
+       -1.01E+02  4.58E+02 -8.63E+00  7.25E+02
 
 TH 5
+        3.65E+01 -1.89E+02 -1.63E+02 -1.06E+02  3.97E+02
 
 TH 6
+        2.81E+00 -1.37E+01  4.36E+00 -2.32E+01 -4.31E-02  2.16E+02
 
 TH 7
+        5.93E-01 -2.32E+00  1.15E-01 -4.26E+00  3.79E+00  8.71E-01  3.26E+00
 
 TH 8
+       -6.09E-03  5.73E-03 -1.67E-02  1.17E-03  2.18E-02  7.32E-02  1.64E-03  3.42E+01
 
 TH 9
+        1.72E+00 -2.38E+01  6.28E+00 -1.25E+01  1.17E+01  2.44E-01  1.26E+01  4.88E-03  1.42E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.01E+01 -1.36E+01 -1.32E+00 -1.80E+01  4.74E+00  4.17E+00  1.48E+00  4.34E-03  1.45E+01  0.00E+00  4.13E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.984
Stop Time:
Wed Sep 29 23:55:04 CDT 2021
