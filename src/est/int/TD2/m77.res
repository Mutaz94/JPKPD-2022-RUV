Sat Sep 25 05:05:46 CDT 2021
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
$DATA ../../../../data/int/TD2/dat77.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3841.80003525030        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.6860E+00  2.8135E+01  1.6427E+01  3.2275E+01 -9.4730E+01 -6.6315E+00 -2.1813E+01 -1.0240E+01  4.0415E+01  2.2847E+00
            -3.8571E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3854.53553085060        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0167E+00  1.0409E+00  1.2605E+00  9.8527E-01  1.1650E+00  1.0042E+00  1.1095E+00  1.3197E+00  7.5864E-01  1.0491E+00
             1.0261E+00
 PARAMETER:  1.1656E-01  1.4013E-01  3.3147E-01  8.5161E-02  2.5273E-01  1.0417E-01  2.0393E-01  3.7737E-01 -1.7622E-01  1.4791E-01
             1.2578E-01
 GRADIENT:   3.7871E+01 -5.5277E+00  2.7351E+01  6.4436E+01 -1.2860E+01 -4.0719E+00  5.7272E+00 -1.8175E+01 -9.4942E+00  7.0646E+00
             1.5920E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3857.48043267337        NO. OF FUNC. EVALS.:  88
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  1.0147E+00  1.0781E+00  1.3521E+00  9.6555E-01  1.2280E+00  9.9920E-01  1.0546E+00  1.4826E+00  7.6588E-01  1.0913E+00
             1.0125E+00
 PARAMETER:  1.1462E-01  1.7518E-01  4.0169E-01  6.4946E-02  3.0538E-01  9.9195E-02  1.5313E-01  4.9380E-01 -1.6674E-01  1.8739E-01
             1.1239E-01
 GRADIENT:   3.4231E+01 -3.9629E+00  2.3481E+01  5.6887E+01 -1.0559E+00 -6.3011E+00  1.1826E+00 -1.5521E+01 -4.3252E+00  5.7955E+00
            -7.8663E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3858.12461699499        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:      286
 NPARAMETR:  1.0135E+00  1.0861E+00  1.3544E+00  9.6018E-01  1.2358E+00  1.0082E+00  1.0455E+00  1.5068E+00  7.6307E-01  1.1092E+00
             1.0129E+00
 PARAMETER:  1.1340E-01  1.8257E-01  4.0338E-01  5.9370E-02  3.1171E-01  1.0815E-01  1.4451E-01  5.0998E-01 -1.7041E-01  2.0366E-01
             1.1278E-01
 GRADIENT:  -2.0390E+01 -1.8293E+01  1.8878E+01  4.4067E+01 -1.4174E+01 -7.9514E+00 -3.0239E-01 -1.4307E+01 -5.5352E+00  6.4266E+00
            -6.6582E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3858.20321715053        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      476
 NPARAMETR:  1.0135E+00  1.0861E+00  1.3545E+00  9.6020E-01  1.2358E+00  1.0283E+00  1.0474E+00  1.5069E+00  7.6307E-01  1.1093E+00
             1.0128E+00
 PARAMETER:  1.1341E-01  1.8259E-01  4.0343E-01  5.9383E-02  3.1175E-01  1.2787E-01  1.4632E-01  5.1004E-01 -1.7041E-01  2.0369E-01
             1.1276E-01
 GRADIENT:  -1.9546E+01 -1.8154E+01  1.8862E+01  4.4087E+01 -1.4155E+01  1.0016E-03 -2.2390E-03 -1.4315E+01 -5.4827E+00  6.5072E+00
            -6.6106E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3858.99201648597        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      636
 NPARAMETR:  1.0225E+00  1.0889E+00  1.3394E+00  9.5661E-01  1.2382E+00  1.0274E+00  1.0472E+00  1.5339E+00  7.6786E-01  1.1017E+00
             1.0161E+00
 PARAMETER:  1.2230E-01  1.8520E-01  3.9223E-01  5.5637E-02  3.1366E-01  1.2706E-01  1.4610E-01  5.2780E-01 -1.6415E-01  1.9684E-01
             1.1594E-01
 GRADIENT:   5.4426E+01 -4.8205E+00  1.4780E+01  4.8341E+01  1.2665E+00  6.4542E+00  1.3941E+00 -9.0217E+00 -3.4969E+00  5.5234E+00
             1.2159E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3859.08480607002        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      766
 NPARAMETR:  1.0224E+00  1.0895E+00  1.3381E+00  9.5605E-01  1.2385E+00  1.0217E+00  1.0468E+00  1.5349E+00  7.9549E-01  1.1013E+00
             1.0161E+00
 PARAMETER:  1.2214E-01  1.8568E-01  3.9125E-01  5.5050E-02  3.1393E-01  1.2151E-01  1.4573E-01  5.2845E-01 -1.2880E-01  1.9646E-01
             1.1596E-01
 GRADIENT:   5.4006E+01 -5.7765E+00  1.5365E+01  4.7655E+01  2.1675E+00  3.9396E+00  2.4355E+00 -8.1225E+00  2.5335E+00  5.7592E+00
             2.8221E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3859.13599548207        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      928            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0224E+00  1.0897E+00  1.3373E+00  9.5580E-01  1.2387E+00  1.0278E+00  1.0467E+00  1.5362E+00  7.8711E-01  1.1008E+00
             1.0161E+00
 PARAMETER:  1.2215E-01  1.8594E-01  3.9064E-01  5.4792E-02  3.1406E-01  1.2738E-01  1.4564E-01  5.2933E-01 -1.3938E-01  1.9606E-01
             1.1594E-01
 GRADIENT:   5.4053E+01 -5.4640E+00  1.4758E+01  4.7375E+01  2.0681E+00  6.5857E+00  2.1288E+00 -8.0942E+00  6.8568E-01  5.5072E+00
             2.3810E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3859.30990707255        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1101
 NPARAMETR:  1.0224E+00  1.0928E+00  1.3373E+00  9.5580E-01  1.2406E+00  1.0273E+00  1.0455E+00  1.5528E+00  7.8699E-01  1.0957E+00
             1.0158E+00
 PARAMETER:  1.2215E-01  1.8870E-01  3.9065E-01  5.4795E-02  3.1557E-01  1.2693E-01  1.4447E-01  5.4004E-01 -1.3954E-01  1.9137E-01
             1.1572E-01
 GRADIENT:  -1.2983E+00 -1.7538E+01  1.0820E+01  4.1151E+01 -1.0689E+01 -2.0265E-01  8.6619E-01 -7.1661E+00  3.5803E-03  3.3359E+00
             1.8357E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -3859.31640790661        NO. OF FUNC. EVALS.:  68
 CUMULATIVE NO. OF FUNC. EVALS.:     1169
 NPARAMETR:  1.0224E+00  1.0929E+00  1.3373E+00  9.5580E-01  1.2407E+00  1.0274E+00  1.0453E+00  1.5533E+00  7.8696E-01  1.0957E+00
             1.0158E+00
 PARAMETER:  1.2215E-01  1.8884E-01  3.9065E-01  5.4795E-02  3.1570E-01  1.2694E-01  1.4432E-01  5.4042E-01 -1.3958E-01  1.9139E-01
             1.1572E-01
 GRADIENT:   6.8488E+05  4.4300E+05  1.0710E+05  4.1833E+05  2.6499E+05 -1.8917E-01  2.8985E+05  1.5466E+05 -3.2386E-03  2.1856E+05
            -7.2299E+05
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.1         3.3         3.3         3.9         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1169
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.4442E-04 -1.4679E-02 -4.3371E-02  2.7231E-03 -2.9084E-02
 SE:             2.9930E-02  2.3677E-02  2.1033E-02  2.5039E-02  2.3907E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7483E-01  5.3527E-01  3.9204E-02  9.1340E-01  2.2378E-01

 ETASHRINKSD(%)  1.0000E-10  2.0680E+01  2.9536E+01  1.6115E+01  1.9909E+01
 ETASHRINKVR(%)  1.0000E-10  3.7083E+01  5.0348E+01  2.9633E+01  3.5854E+01
 EBVSHRINKSD(%)  2.2626E-01  2.0818E+01  3.4767E+01  1.8278E+01  1.6838E+01
 EBVSHRINKVR(%)  4.5201E-01  3.7302E+01  5.7446E+01  3.3215E+01  3.0841E+01
 RELATIVEINF(%)  9.9546E+01  3.1169E+01  3.4209E+01  3.5856E+01  3.6310E+01
 EPSSHRINKSD(%)  2.1460E+01
 EPSSHRINKVR(%)  3.8315E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3859.3164079066123     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2205.2270481382016     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3859.316       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.09E+00  1.34E+00  9.56E-01  1.24E+00  1.03E+00  1.05E+00  1.55E+00  7.87E-01  1.10E+00  1.02E+00
 


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
+        1.34E+09
 
 TH 2
+       -2.74E+04  4.91E+08
 
 TH 3
+        2.65E+02 -6.51E+03  7.67E+07
 
 TH 4
+        1.44E+03 -3.53E+04 -1.67E+04  2.29E+09
 
 TH 5
+        3.53E+02 -8.88E+03 -3.01E+03  5.59E+08  1.36E+08
 
 TH 6
+        1.40E+03  8.48E+02  3.34E+02  1.83E+03  4.44E+02  1.85E+02
 
 TH 7
+        1.11E+09  6.72E+08 -3.29E+03 -1.80E+04 -4.38E+03  1.16E+03  9.19E+08
 
 TH 8
+       -3.49E+05  1.21E+08 -8.35E+04 -4.56E+05 -1.11E+05  2.09E+02  1.65E+08  2.96E+07
 
 TH 9
+       -3.68E+04 -2.23E+04 -8.78E+03 -4.80E+04 -1.17E+04  4.11E-01 -3.04E+04 -5.45E+03  1.42E+02
 
 TH10
+       -6.22E+03  4.83E+08 -1.48E+03 -8.10E+03 -2.01E+03  8.34E+02  6.61E+08  1.19E+08 -2.19E+04  4.76E+08
 
 TH11
+       -1.42E+09  2.90E+04 -1.72E+04 -1.86E+09 -2.29E+04 -1.49E+03  1.46E+04  3.71E+05  3.91E+04  6.61E+03  1.51E+09
 
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
 #CPUT: Total CPU Time in Seconds,       49.122
Stop Time:
Sat Sep 25 05:06:36 CDT 2021
