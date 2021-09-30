Wed Sep 29 09:35:16 CDT 2021
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
$DATA ../../../../data/int/D/dat71.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26168.6079223803        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6205E+02  3.6655E+02 -2.8891E+01 -2.9521E+01  2.7845E+02 -2.6091E+03 -1.2360E+03 -9.1968E+01 -2.2792E+03 -7.9628E+02
            -5.2590E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -927.436378394527        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  2.3733E+00  1.7188E+00  7.6743E-01  3.5272E+00  9.4912E-01  6.8838E+00  6.1963E+00  1.0602E+00  7.0746E+00  2.8461E+00
             9.4570E+00
 PARAMETER:  9.6427E-01  6.4162E-01 -1.6470E-01  1.3605E+00  4.7782E-02  2.0292E+00  1.9240E+00  1.5848E-01  2.0565E+00  1.1459E+00
             2.3468E+00
 GRADIENT:   8.3378E+01  9.3542E+00 -3.5270E+01  6.8183E+01 -6.8990E+01  3.6162E+02  2.0070E+02  5.3771E+00  1.4914E+02  5.9793E+01
             3.1726E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1019.91090200840        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  1.0801E+00  4.9672E+00  2.4255E+01  2.9243E+00  3.4945E+00  9.0519E+00  9.1998E+00  1.1250E+00  1.0608E+01  1.9260E+00
             8.9579E+00
 PARAMETER:  1.7706E-01  1.7029E+00  3.2886E+00  1.1731E+00  1.3512E+00  2.3030E+00  2.3192E+00  2.1776E-01  2.4617E+00  7.5545E-01
             2.2925E+00
 GRADIENT:  -4.3154E+00  1.6150E+01  3.6285E+00  1.9805E+01  4.5181E+01  1.8147E+02  8.9678E+01 -1.1181E+01  3.7423E+01  3.5056E+01
             2.4291E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1174.44185873874        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.7397E+00  1.4379E+00  2.0396E+01  2.2528E+00  2.2541E+00  3.8759E+00  7.1396E+00  8.6192E+00  5.0440E+00  1.1772E+00
             9.6956E+00
 PARAMETER:  6.5373E-01  4.6317E-01  3.1153E+00  9.1220E-01  9.1275E-01  1.4548E+00  2.0657E+00  2.2540E+00  1.7182E+00  2.6316E-01
             2.3717E+00
 GRADIENT:   4.1577E+01 -4.6091E-01 -1.8547E+01  2.9952E+01 -6.8622E+01  5.4836E+01  2.8638E+00  5.1719E+01  6.6990E+01  1.4820E+01
             3.2453E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1270.48015113570        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  9.5398E-01  2.4010E+00  1.5018E+01  5.7685E-01  2.6713E+00  2.9191E+00  4.6280E+00  5.3945E+00  3.0541E+00  8.2606E-01
             7.7514E+00
 PARAMETER:  5.2883E-02  9.7587E-01  2.8093E+00 -4.5018E-01  1.0826E+00  1.1713E+00  1.6321E+00  1.7854E+00  1.2165E+00 -9.1093E-02
             2.1479E+00
 GRADIENT:  -5.0755E+01 -5.6568E+00  1.1515E+00 -1.2965E+01 -1.2748E-03 -2.5338E+01  9.7053E+00 -5.3510E-01  1.0631E+01  3.5992E+00
             6.4179E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1285.82706218809        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      809             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3021E+00  1.7233E+00  1.3417E+01  8.8844E-01  2.5866E+00  3.2297E+00  5.3574E+00  4.6839E+00  3.0539E+00  7.0325E-01
             7.6036E+00
 PARAMETER:  3.6397E-01  6.4424E-01  2.6965E+00 -1.8292E-02  1.0503E+00  1.2724E+00  1.7785E+00  1.6441E+00  1.2164E+00 -2.5204E-01
             2.1286E+00
 GRADIENT:   4.3566E+01  1.2056E+01 -2.0519E-01 -1.0446E+01  1.4030E+01  1.4612E+02  1.5362E+02  1.2433E+01  3.0188E+01  3.6197E+00
             2.0189E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1290.33097994920        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      927            RESET HESSIAN, TYPE II
 NPARAMETR:  1.2916E+00  1.7497E+00  1.3702E+01  9.0058E-01  2.5641E+00  3.1960E+00  5.0326E+00  3.6069E+00  2.8498E+00  6.9694E-01
             7.4707E+00
 PARAMETER:  3.5586E-01  6.5942E-01  2.7175E+00 -4.7122E-03  1.0416E+00  1.2619E+00  1.7159E+00  1.3828E+00  1.1472E+00 -2.6105E-01
             2.1110E+00
 GRADIENT:   4.3130E+01  1.2531E+01 -1.8769E-01 -4.7263E+00  7.4238E+00  1.4613E+02  1.3343E+02  2.5934E+00  2.2705E+01  1.5507E+00
            -1.5660E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1291.50006634562        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  1.2818E+00  1.7744E+00  1.3736E+01  9.1154E-01  2.5622E+00  3.2024E+00  5.0277E+00  3.5993E+00  2.7072E+00  6.8601E-01
             7.4870E+00
 PARAMETER:  3.4828E-01  6.7345E-01  2.7200E+00  7.3765E-03  1.0409E+00  1.2639E+00  1.7150E+00  1.3807E+00  1.0959E+00 -2.7686E-01
             2.1132E+00
 GRADIENT:   1.1670E+01 -5.3583E+00  1.5768E-01 -5.1807E+00  2.3418E+00  1.3229E+00 -2.8647E+01  1.4291E+00  1.3461E+01  1.4544E+00
            -5.2151E+01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1291.51483794395        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1137
 NPARAMETR:  1.2838E+00  1.7706E+00  1.3904E+01  9.1204E-01  2.5530E+00  3.2222E+00  4.9907E+00  3.5492E+00  2.7202E+00  6.8365E-01
             7.4340E+00
 PARAMETER:  3.4828E-01  6.7433E-01  2.7200E+00  8.3765E-03  1.0408E+00  1.2639E+00  1.7150E+00  1.3805E+00  1.0959E+00 -2.7905E-01
             2.1132E+00
 GRADIENT:  -1.7706E+03  9.1716E+02 -2.2873E+02  6.2159E+03  2.1395E+00 -2.5940E+02  3.4818E+02  1.3381E+00 -5.5429E+02  1.1158E+03
             2.2427E+02
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.4         2.3         2.3         1.9         2.3         2.3
                    2.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1137
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.0965E-02  3.2442E-02 -9.1751E-03 -5.1018E-02 -8.0200E-03
 SE:             3.1100E-02  2.6896E-02  4.0293E-03  1.1496E-02  9.8445E-03
 N:                     100         100         100         100         100

 P VAL.:         7.2440E-01  2.2773E-01  2.2781E-02  9.0902E-06  4.1527E-01

 ETASHRINKSD(%)  1.0000E-10  9.8957E+00  8.6501E+01  6.1487E+01  6.7020E+01
 ETASHRINKVR(%)  1.0000E-10  1.8812E+01  9.8178E+01  8.5168E+01  8.9123E+01
 EBVSHRINKSD(%)  1.5978E+00  1.0577E+01  8.9208E+01  5.6137E+01  6.6304E+01
 EBVSHRINKVR(%)  3.1701E+00  2.0035E+01  9.8835E+01  8.0761E+01  8.8646E+01
 RELATIVEINF(%)  9.6795E+01  2.9687E+01  4.5224E-01  7.0231E+00  4.6643E+00
 EPSSHRINKSD(%)  7.6514E+00
 EPSSHRINKVR(%)  1.4717E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1291.5148379439488     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       362.57452182446195     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.27
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    19.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1291.515       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.28E+00  1.78E+00  1.37E+01  9.12E-01  2.56E+00  3.20E+00  5.03E+00  3.60E+00  2.71E+00  6.85E-01  7.49E+00
 


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
+        7.77E+04
 
 TH 2
+        1.11E+01  1.09E+04
 
 TH 3
+       -6.84E-02  1.89E-01  1.13E+01
 
 TH 4
+       -5.85E-01  1.15E+01 -9.67E-01  1.87E+06
 
 TH 5
+        1.53E+02 -5.92E+01  2.40E+00 -7.57E+02  4.32E+01
 
 TH 6
+       -4.45E+01  1.25E+00  1.57E-02  4.89E-02  1.70E+01  9.01E+02
 
 TH 7
+        1.93E+01  2.02E+00 -4.10E-02 -7.10E+00 -7.47E+00  3.58E+00  2.13E+02
 
 TH 8
+        3.71E-01 -1.18E+00 -9.74E-01  7.15E+00 -8.89E+00 -1.63E-01  4.98E-01  8.24E+00
 
 TH 9
+       -5.58E+01  2.32E-01 -2.02E-01 -5.43E-01  2.22E+01 -1.62E-01  1.08E+00  1.29E+00  1.77E+03
 
 TH10
+        1.04E+02 -3.88E+01 -2.18E+03 -5.20E+02  3.04E+04  1.18E+01 -5.37E+00 -7.35E+00  1.50E+01  4.26E+05
 
 TH11
+        7.59E+00 -1.25E+00  2.81E-01 -5.19E+00 -2.87E+00  3.10E+01 -8.80E-01 -2.48E+00  3.60E-01  1.39E+00  7.41E+01
 
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
 #CPUT: Total CPU Time in Seconds,       63.781
Stop Time:
Wed Sep 29 09:36:21 CDT 2021
