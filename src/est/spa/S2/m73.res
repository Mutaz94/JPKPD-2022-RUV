Sat Sep 18 13:36:59 CDT 2021
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
$DATA ../../../../data/spa/S2/dat73.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m73.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1693.37007940631        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.6502E+00 -1.5377E+01  4.1546E+01 -7.8027E+01 -8.2807E+01  5.3595E+01 -8.1684E+00 -1.4828E+01 -1.6842E+01  3.5559E+01
            -1.9607E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1709.00543268939        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0461E+00  1.1192E+00  9.9226E-01  9.2522E-01  1.1485E+00  8.1194E-01  1.1607E+00  1.2895E+00  1.0953E+00  6.0959E-01
             1.0131E+00
 PARAMETER:  1.4512E-01  2.1262E-01  9.2234E-02  2.2280E-02  2.3848E-01 -1.0833E-01  2.4899E-01  3.5428E-01  1.9105E-01 -3.9497E-01
             1.1302E-01
 GRADIENT:   1.5646E+02 -5.5527E+01 -2.1301E+01 -2.7838E+01  8.7656E+01 -2.4599E+01  2.2885E+00 -1.0289E+01  4.5173E+00 -4.1919E+00
            -4.0726E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1710.25298027965        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  1.0557E+00  1.2467E+00  9.9386E-01  9.1282E-01  1.0691E+00  8.6831E-01  1.1269E+00  1.7374E+00  1.0391E+00  4.0882E-01
             1.0101E+00
 PARAMETER:  1.5418E-01  3.2052E-01  9.3837E-02  8.7853E-03  1.6682E-01 -4.1203E-02  2.1945E-01  6.5241E-01  1.3837E-01 -7.9449E-01
             1.1007E-01
 GRADIENT:   1.7380E+02  6.8430E+01  3.0152E+01  1.1767E+01 -8.1902E+01  1.4948E+00  6.1651E+00  2.9049E+00 -4.3032E+00 -2.4803E+00
            -6.0440E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1715.55214157524        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  1.0047E+00  1.1650E+00  1.0693E+00  9.4141E-01  1.1234E+00  8.5161E-01  1.0809E+00  1.6620E+00  1.0739E+00  5.6182E-01
             1.0140E+00
 PARAMETER:  1.0468E-01  2.5272E-01  1.6697E-01  3.9624E-02  2.1639E-01 -6.0624E-02  1.7778E-01  6.0800E-01  1.7132E-01 -4.7656E-01
             1.1394E-01
 GRADIENT:   8.6351E+00  9.6950E+00  4.0138E+00 -5.6036E-01 -9.3164E+00 -1.8865E+00  1.9454E+00 -1.5196E-01 -1.3003E-01  1.4776E-02
             7.4042E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1715.84592722600        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      428
 NPARAMETR:  1.0092E+00  1.1853E+00  1.1147E+00  9.3321E-01  1.1540E+00  8.5575E-01  1.0522E+00  1.7705E+00  1.0985E+00  5.8623E-01
             1.0157E+00
 PARAMETER:  1.0916E-01  2.6997E-01  2.0856E-01  3.0870E-02  2.4325E-01 -5.5778E-02  1.5084E-01  6.7127E-01  1.9392E-01 -4.3405E-01
             1.1562E-01
 GRADIENT:  -2.1182E+01 -1.2425E+00  2.0742E+00 -2.9323E+00 -5.2847E+00 -2.9004E+00  2.3353E+00  1.6917E-01 -2.4932E-01 -1.4923E-01
             1.1782E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1715.99008087976        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:      575
 NPARAMETR:  1.0160E+00  1.1852E+00  1.1147E+00  9.3511E-01  1.1570E+00  8.6177E-01  1.0195E+00  1.7702E+00  1.1112E+00  5.8616E-01
             1.0141E+00
 PARAMETER:  1.1590E-01  2.6990E-01  2.0861E-01  3.2912E-02  2.4586E-01 -4.8768E-02  1.1932E-01  6.7110E-01  2.0547E-01 -4.3416E-01
             1.1396E-01
 GRADIENT:  -1.0011E+00 -2.2075E+00  1.6019E-01  3.4204E-01 -3.0462E-01 -3.2906E-02  6.2812E-01 -4.0102E-01 -2.4973E-01 -8.7198E-01
             1.0027E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1716.00646067262        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  1.0163E+00  1.1863E+00  1.1148E+00  9.3450E-01  1.1573E+00  8.6184E-01  1.0034E+00  1.7740E+00  1.1195E+00  5.9034E-01
             1.0138E+00
 PARAMETER:  1.1619E-01  2.7086E-01  2.0868E-01  3.2260E-02  2.4612E-01 -4.8685E-02  1.0337E-01  6.7324E-01  2.1285E-01 -4.2705E-01
             1.1373E-01
 GRADIENT:   4.8385E+01  1.0391E+01  2.8396E-01  6.2384E+00 -3.6752E-01  3.1300E+00  5.3319E-01 -7.2813E-02  1.2310E+00 -6.2874E-01
             3.4148E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1716.01319300857        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      899
 NPARAMETR:  1.0163E+00  1.1875E+00  1.1148E+00  9.3421E-01  1.1579E+00  8.6183E-01  1.0027E+00  1.7743E+00  1.1195E+00  5.9376E-01
             1.0136E+00
 PARAMETER:  1.1618E-01  2.7182E-01  2.0867E-01  3.1944E-02  2.4657E-01 -4.8703E-02  1.0270E-01  6.7338E-01  2.1292E-01 -4.2128E-01
             1.1346E-01
 GRADIENT:  -2.1913E-01 -1.8498E+00  4.2583E-02  1.0813E+00 -2.0568E+00 -9.1569E-03  1.6672E-01 -3.7871E-01 -9.4194E-02 -4.8387E-01
             2.6656E-01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1716.01319300857        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:      927
 NPARAMETR:  1.0163E+00  1.1875E+00  1.1148E+00  9.3421E-01  1.1579E+00  8.6183E-01  1.0027E+00  1.7743E+00  1.1195E+00  5.9376E-01
             1.0136E+00
 PARAMETER:  1.1618E-01  2.7182E-01  2.0867E-01  3.1944E-02  2.4657E-01 -4.8703E-02  1.0270E-01  6.7338E-01  2.1292E-01 -4.2128E-01
             1.1346E-01
 GRADIENT:  -1.5323E-01  2.0364E+05 -2.6524E+05 -5.5355E+05  2.2449E+05 -1.0557E-02  2.6949E+05  8.2107E+04  2.5997E+05  1.3139E+05
            -4.8788E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      927
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0835E-03 -1.6338E-02 -3.4473E-02  1.0926E-02 -4.1611E-02
 SE:             2.9834E-02  1.9403E-02  1.8810E-02  2.4458E-02  1.5233E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7103E-01  3.9976E-01  6.6848E-02  6.5508E-01  6.3004E-03

 ETASHRINKSD(%)  5.2274E-02  3.4997E+01  3.6985E+01  1.8061E+01  4.8969E+01
 ETASHRINKVR(%)  1.0452E-01  5.7746E+01  6.0291E+01  3.2860E+01  7.3958E+01
 EBVSHRINKSD(%)  5.8823E-01  3.4866E+01  3.6447E+01  1.9131E+01  4.9694E+01
 EBVSHRINKVR(%)  1.1730E+00  5.7575E+01  5.9611E+01  3.4603E+01  7.4693E+01
 RELATIVEINF(%)  9.7931E+01  1.2424E+00  3.2966E+00  2.2924E+00  5.9892E+00
 EPSSHRINKSD(%)  4.3906E+01
 EPSSHRINKVR(%)  6.8534E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1716.0131930085679     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -980.86236644482972     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.39
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1716.013       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.19E+00  1.11E+00  9.34E-01  1.16E+00  8.62E-01  1.00E+00  1.77E+00  1.12E+00  5.94E-01  1.01E+00
 


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
+        9.93E+08
 
 TH 2
+        1.40E+03  1.33E+08
 
 TH 3
+       -1.96E+03  1.20E+04  2.56E+08
 
 TH 4
+       -4.90E+03 -2.29E+04  3.23E+04  1.59E+09
 
 TH 5
+       -4.11E+08  5.07E+03 -7.51E+03 -5.19E+08  1.70E+08
 
 TH 6
+        1.63E-01  1.10E+03 -1.53E+03 -3.81E+03  1.25E+03  2.67E+02
 
 TH 7
+       -1.14E+09  3.03E+03 -4.19E+03 -1.44E+09  4.71E+08  3.46E+03  1.30E+09
 
 TH 8
+        3.83E+02 -2.32E+03  7.06E+03 -6.28E+03  1.41E+03  2.98E+02  8.23E+02  9.67E+06
 
 TH 9
+       -4.92E+08 -8.85E+03  1.23E+04 -6.21E+08  2.03E+08  1.49E+03  5.64E+08 -2.38E+03  2.44E+08
 
 TH10
+        1.83E+03 -1.10E+04  1.65E+03 -3.01E+04  6.70E+03  1.42E+03  3.92E+03 -3.11E+02 -1.14E+04  2.21E+08
 
 TH11
+       -3.99E+03 -1.41E+04  1.96E+04  6.54E+04 -1.48E+04 -3.09E+03 -8.49E+03 -3.80E+03  2.48E+04 -1.82E+04  1.05E+09
 
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
 #CPUT: Total CPU Time in Seconds,       17.633
Stop Time:
Sat Sep 18 13:37:18 CDT 2021
