Wed Sep 29 13:44:25 CDT 2021
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
$DATA ../../../../data/spa/A3/dat69.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   542.447253187372        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1058E+02  1.3243E+02  1.0702E+02  6.7489E+01  2.0578E+02  3.7841E+01 -5.6674E+01 -4.2176E+01 -1.1612E+02 -2.1475E+02
            -3.9493E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1201.37850841872        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0917E+00  9.5872E-01  9.1262E-01  1.1244E+00  9.3263E-01  8.4360E-01  9.5300E-01  9.8678E-01  9.9189E-01  1.0478E+00
             5.3136E+00
 PARAMETER:  1.8774E-01  5.7841E-02  8.5655E-03  2.1726E-01  3.0256E-02 -7.0072E-02  5.1860E-02  8.6689E-02  9.1855E-02  1.4673E-01
             1.7703E+00
 GRADIENT:   7.5149E+01 -8.6415E+00 -1.7586E+01  8.4746E+00  2.6588E+00 -1.3478E+01  1.0632E+01  6.1017E+00  1.9685E+01  2.2419E+01
             1.4228E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1215.12586103120        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0680E+00  7.2893E-01  3.0903E-01  1.1344E+00  3.9839E-01  9.1339E-01  3.0956E-01  2.1150E-01  1.1279E+00  5.2322E-01
             4.8434E+00
 PARAMETER:  1.6581E-01 -2.1617E-01 -1.0743E+00  2.2607E-01 -8.2033E-01  9.4114E-03 -1.0726E+00 -1.4535E+00  2.2034E-01 -5.4776E-01
             1.6776E+00
 GRADIENT:   2.1217E+00  7.8946E+01  3.1683E+01  8.2651E+01 -8.1828E+01 -6.4102E+00  8.1834E-02  4.8152E-01  1.2549E+01  7.8780E+00
             1.0390E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1223.41956893781        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0233E+00  5.6213E-01  1.7215E-01  1.1103E+00  2.8211E-01  9.7150E-01  1.8658E-01  3.4032E-02  1.1859E+00  4.6167E-01
             4.0584E+00
 PARAMETER:  1.2306E-01 -4.7601E-01 -1.6594E+00  2.0463E-01 -1.1655E+00  7.1084E-02 -1.5789E+00 -3.2804E+00  2.7048E-01 -6.7290E-01
             1.5008E+00
 GRADIENT:  -4.0711E+01  6.2494E+01 -8.3546E+00  1.6001E+02 -1.8522E+01 -3.9681E+00 -2.8081E-01  1.4127E-04 -2.0999E+01  7.4359E-01
             1.7428E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1233.44128178909        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      406
 NPARAMETR:  1.0534E+00  5.3728E-01  2.6525E-01  1.1060E+00  3.4394E-01  9.4188E-01  3.3371E-01  4.3205E-02  1.0708E+00  3.1305E-01
             4.0774E+00
 PARAMETER:  1.5205E-01 -5.2123E-01 -1.2271E+00  2.0079E-01 -9.6728E-01  4.0118E-02 -9.9749E-01 -3.0418E+00  1.6844E-01 -1.0614E+00
             1.5055E+00
 GRADIENT:   1.6462E+00  4.6471E+00 -7.0922E+00  2.3335E+01  1.0823E+01 -3.4150E+00 -8.5331E-01  7.2807E-03 -1.4835E+00 -1.5118E+00
            -1.5593E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1234.66545428431        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      581
 NPARAMETR:  1.0531E+00  4.7206E-01  2.3965E-01  1.0887E+00  3.0805E-01  9.5686E-01  6.7377E-01  1.9529E-02  1.0683E+00  2.5367E-01
             4.0761E+00
 PARAMETER:  1.5169E-01 -6.5064E-01 -1.3286E+00  1.8496E-01 -1.0775E+00  5.5907E-02 -2.9486E-01 -3.8358E+00  1.6606E-01 -1.2717E+00
             1.5051E+00
 GRADIENT:   7.1057E-02  2.6986E-02 -1.1308E+00  5.7412E+00  2.9949E+00  1.5002E-01 -3.3563E-01 -2.8186E-04 -1.9129E+00 -1.2842E+00
            -5.5092E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1235.49503534997        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      761
 NPARAMETR:  1.0415E+00  5.4308E-01  1.7550E-01  9.8534E-01  2.8380E-01  9.6848E-01  5.9871E-01  1.5219E-02  1.1760E+00  3.0298E-01
             3.9785E+00
 PARAMETER:  1.4068E-01 -5.1050E-01 -1.6401E+00  8.5228E-02 -1.1595E+00  6.7977E-02 -4.1298E-01 -4.0852E+00  2.6216E-01 -1.0941E+00
             1.4809E+00
 GRADIENT:  -4.7690E+00  5.8745E+00  1.0752E+00  8.8047E+00 -7.4830E+00 -3.5592E-01 -1.3640E-01 -1.1073E-03 -2.6086E+00 -1.2053E+00
            -2.4777E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1235.74921136776        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      936
 NPARAMETR:  1.0401E+00  5.3316E-01  1.6278E-01  9.6599E-01  2.7365E-01  9.7262E-01  4.9540E-01  1.0000E-02  1.2242E+00  3.9700E-01
             3.9232E+00
 PARAMETER:  1.3933E-01 -5.2894E-01 -1.7153E+00  6.5393E-02 -1.1959E+00  7.2236E-02 -6.0238E-01 -4.7268E+00  3.0227E-01 -8.2382E-01
             1.4669E+00
 GRADIENT:  -1.6228E+00 -3.4837E+00 -5.6407E+00  5.1930E+00  9.8442E+00 -3.7447E-01  2.1242E-01  0.0000E+00 -8.4244E-01  3.4635E-01
            -1.3901E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1235.80089078295        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  1.0427E+00  5.2171E-01  1.6579E-01  9.6968E-01  2.7140E-01  9.7314E-01  4.7786E-01  1.0000E-02  1.2220E+00  4.0213E-01
             3.9364E+00
 PARAMETER:  1.4183E-01 -5.5065E-01 -1.6970E+00  6.9212E-02 -1.2042E+00  7.2768E-02 -6.3845E-01 -4.7681E+00  3.0052E-01 -8.1097E-01
             1.4703E+00
 GRADIENT:   7.2862E-02  5.0145E-01  1.7527E-01  4.7050E-01 -5.0751E-01  2.2138E-03  2.6098E-02  0.0000E+00 -1.1290E-02  6.8874E-02
             2.4008E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1028
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.2414E-04 -2.9405E-03  8.1618E-05 -1.5263E-02  6.6438E-03
 SE:             2.8374E-02  8.8134E-03  1.5194E-04  2.4781E-02  1.3599E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9370E-01  7.3865E-01  5.9114E-01  5.3794E-01  6.2515E-01

 ETASHRINKSD(%)  4.9419E+00  7.0474E+01  9.9491E+01  1.6982E+01  5.4443E+01
 ETASHRINKVR(%)  9.6397E+00  9.1282E+01  9.9997E+01  3.1080E+01  7.9245E+01
 EBVSHRINKSD(%)  4.6403E+00  7.0587E+01  9.9500E+01  1.5289E+01  5.4370E+01
 EBVSHRINKVR(%)  9.0653E+00  9.1349E+01  9.9997E+01  2.8240E+01  7.9179E+01
 RELATIVEINF(%)  8.1501E+01  2.6216E-01  2.0156E-04  2.4927E+01  3.8972E-01
 EPSSHRINKSD(%)  2.4688E+01
 EPSSHRINKVR(%)  4.3281E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1235.8008907829544     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -500.65006421921623     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1235.801       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  5.22E-01  1.66E-01  9.70E-01  2.71E-01  9.73E-01  4.78E-01  1.00E-02  1.22E+00  4.02E-01  3.94E+00
 


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
+       -9.19E+01  1.48E+03
 
 TH 3
+       -3.94E+02  2.95E+03  9.75E+03
 
 TH 4
+       -5.94E+01  1.73E+02 -8.40E+02  5.66E+02
 
 TH 5
+        5.05E+02 -5.10E+03 -1.19E+04 -1.02E+02  1.92E+04
 
 TH 6
+       -2.56E+00 -4.79E+00  5.86E+01 -1.76E+01  3.80E+01  1.72E+02
 
 TH 7
+       -1.64E+00 -2.96E+01 -7.47E+01 -1.91E+00  1.14E+02  9.36E-01  9.33E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.76E+00 -2.51E+01  1.10E+02 -1.06E+01  7.94E+01  1.44E+00  3.06E+00  0.00E+00  6.22E+01
 
 TH10
+       -8.69E+00 -8.90E+01 -2.33E+02 -3.64E+00  3.96E+02  2.88E+00  1.64E+01  0.00E+00  1.18E+00  6.44E+01
 
 TH11
+       -1.65E+01 -1.04E+01 -2.67E+01 -6.16E+00  2.84E+01  2.89E+00  4.95E+00  0.00E+00  6.41E+00  1.59E+01  2.37E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.997
Stop Time:
Wed Sep 29 13:44:47 CDT 2021
