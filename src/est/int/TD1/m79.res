Sat Sep 25 04:21:47 CDT 2021
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
$DATA ../../../../data/int/TD1/dat79.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3879.23399615303        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3614E+00  2.5445E+01 -2.7018E+01 -5.6257E+01 -1.0035E+01  2.6409E+01  1.2999E+01 -1.0687E+01 -4.6280E+00  3.7866E+00
            -5.0004E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3882.97887070673        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0322E+00  9.5968E-01  1.0382E+00  1.0520E+00  9.8808E-01  8.7659E-01  9.0416E-01  1.1079E+00  9.8420E-01  9.6160E-01
             1.0360E+00
 PARAMETER:  1.3167E-01  5.8849E-02  1.3746E-01  1.5068E-01  8.8010E-02 -3.1716E-02 -7.5246E-04  2.0250E-01  8.4073E-02  6.0840E-02
             1.3540E-01
 GRADIENT:   8.7920E+01 -1.7203E-01  4.6046E-01  1.4863E+01 -5.0233E+00 -2.5529E+01 -1.4210E+00 -1.2985E+00 -1.3606E+00 -5.3620E+00
             2.9316E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3883.24841768136        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      166
 NPARAMETR:  1.0312E+00  9.6494E-01  1.0456E+00  1.0434E+00  1.0016E+00  8.8129E-01  8.6866E-01  1.1258E+00  9.8068E-01  1.0001E+00
             1.0386E+00
 PARAMETER:  1.3072E-01  6.4314E-02  1.4456E-01  1.4250E-01  1.0158E-01 -2.6371E-02 -4.0800E-02  2.1845E-01  8.0496E-02  1.0015E-01
             1.3786E-01
 GRADIENT:   8.4178E+01 -1.2088E+01 -3.9270E-01  3.7212E+00  4.3854E+00 -2.3021E+01 -2.2503E+00 -1.1042E-01 -3.7953E+00 -1.6016E+00
             3.3700E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3884.39539093839        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      239
 NPARAMETR:  1.0124E+00  9.8336E-01  1.0551E+00  1.0340E+00  1.0153E+00  9.1417E-01  8.5716E-01  1.1190E+00  1.0026E+00  1.0056E+00
             1.0277E+00
 PARAMETER:  1.1233E-01  8.3225E-02  1.5364E-01  1.3344E-01  1.1518E-01  1.0259E-02 -5.4133E-02  2.1239E-01  1.0256E-01  1.0554E-01
             1.2733E-01
 GRADIENT:   2.7112E+01 -1.2858E+00  1.7766E+00  1.3767E+00  4.6005E+00 -7.1226E+00 -2.4116E+00 -2.2111E+00  1.0707E+00 -3.2656E+00
             1.0313E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3884.52078437579        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  1.0136E+00  9.8489E-01  1.0569E+00  1.0344E+00  1.0168E+00  9.1821E-01  8.6114E-01  1.1267E+00  1.0021E+00  1.0076E+00
             1.0267E+00
 PARAMETER:  1.1352E-01  8.4770E-02  1.5538E-01  1.3381E-01  1.1670E-01  1.4674E-02 -4.9495E-02  2.1933E-01  1.0213E-01  1.0758E-01
             1.2635E-01
 GRADIENT:  -1.8140E+01 -6.8788E+00  5.7082E-01 -9.4698E+00 -8.2706E-01 -8.7579E+00 -2.3634E+00 -1.9654E+00  6.6833E-02 -3.0464E+00
             8.3654E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3884.65338564985        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      568
 NPARAMETR:  1.0136E+00  9.8496E-01  1.0570E+00  1.0344E+00  1.0169E+00  9.3830E-01  8.7712E-01  1.1269E+00  1.0012E+00  1.0176E+00
             1.0267E+00
 PARAMETER:  1.1353E-01  8.4842E-02  1.5548E-01  1.3379E-01  1.1672E-01  3.6314E-02 -3.1113E-02  2.1948E-01  1.0118E-01  1.1740E-01
             1.2635E-01
 GRADIENT:  -1.7251E+01 -6.5385E+00  2.8763E-01 -9.4181E+00 -1.5274E+00 -6.4074E-02 -1.3419E-01 -1.9248E+00  3.8466E-02  1.8456E-01
             9.0710E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3884.71901502966        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.0150E+00  9.8624E-01  1.0570E+00  1.0358E+00  1.0171E+00  9.3794E-01  8.7993E-01  1.1361E+00  1.0006E+00  1.0151E+00
             1.0249E+00
 PARAMETER:  1.1490E-01  8.6144E-02  1.5539E-01  1.3522E-01  1.1691E-01  3.5930E-02 -2.7913E-02  2.2763E-01  1.0056E-01  1.1502E-01
             1.2463E-01
 GRADIENT:   3.6149E+01  3.5027E+00  5.5951E-01  6.9950E+00  3.6443E+00  3.2957E+00  4.4128E-01 -1.1058E+00  9.3120E-01  4.1095E-01
             6.1552E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3884.72026481892        NO. OF FUNC. EVALS.: 210
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  1.0152E+00  9.8631E-01  1.0571E+00  1.0359E+00  1.0171E+00  9.3837E-01  8.7991E-01  1.1366E+00  1.0005E+00  1.0151E+00
             1.0249E+00
 PARAMETER:  1.1496E-01  8.6173E-02  1.5547E-01  1.3520E-01  1.1691E-01  3.6384E-02 -2.7936E-02  2.2780E-01  1.0053E-01  1.1500E-01
             1.2461E-01
 GRADIENT:  -1.3593E+01 -3.4943E+00 -5.0679E-01 -5.9385E+00 -1.9714E+00 -8.8630E-03 -3.1346E-02 -1.2821E+00 -2.8416E-02 -8.1236E-02
             5.8370E+00
 NUMSIGDIG:         1.0         1.3         2.0         1.3         1.7         3.7         2.3         1.0         2.7         2.2
                    1.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      915
 NO. OF SIG. DIGITS IN FINAL EST.:  1.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2065E-03 -2.8901E-02 -2.7459E-02  1.5606E-02 -2.9312E-02
 SE:             2.9875E-02  2.0979E-02  1.9599E-02  2.7870E-02  2.4853E-02
 N:                     100         100         100         100         100

 P VAL.:         8.3543E-01  1.6833E-01  1.6120E-01  5.7550E-01  2.3823E-01

 ETASHRINKSD(%)  1.0000E-10  2.9717E+01  3.4341E+01  6.6314E+00  1.6739E+01
 ETASHRINKVR(%)  1.0000E-10  5.0603E+01  5.6889E+01  1.2823E+01  3.0677E+01
 EBVSHRINKSD(%)  3.1175E-01  2.9394E+01  3.6642E+01  6.9984E+00  1.6905E+01
 EBVSHRINKVR(%)  6.2253E-01  5.0149E+01  5.9857E+01  1.3507E+01  3.0952E+01
 RELATIVEINF(%)  9.9375E+01  2.2783E+01  2.6755E+01  5.9349E+01  2.7819E+01
 EPSSHRINKSD(%)  2.1745E+01
 EPSSHRINKVR(%)  3.8762E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3884.7202648189186     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2230.6309050505079     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.92
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
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
 





 #OBJV:********************************************    -3884.720       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  9.86E-01  1.06E+00  1.04E+00  1.02E+00  9.38E-01  8.80E-01  1.14E+00  1.00E+00  1.02E+00  1.02E+00
 


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
+        1.22E+03
 
 TH 2
+       -1.14E+00  7.37E+02
 
 TH 3
+        9.26E-01 -1.34E+01  3.03E+02
 
 TH 4
+       -1.10E+00 -2.16E+09 -1.30E+09  1.52E+09
 
 TH 5
+       -3.98E-01 -3.79E+02 -2.11E+02  1.80E+09  7.26E+02
 
 TH 6
+        1.85E+01  2.37E+00  5.28E-01 -1.44E+00  2.33E+00  2.24E+02
 
 TH 7
+       -1.75E+00  1.90E+01 -1.84E+01  2.62E+00 -8.29E+00  1.46E+00  6.32E+01
 
 TH 8
+       -2.02E-01 -1.17E+09 -7.02E+08  8.24E+08  9.71E+08  1.32E+00 -9.11E-01  4.46E+08
 
 TH 9
+        2.56E+00 -2.20E+01  1.17E+01  3.17E+01 -1.29E+00  3.36E+00  7.07E+00  3.37E+00  1.51E+02
 
 TH10
+       -3.03E+00 -2.61E+01  1.03E+01  3.92E+00 -2.10E+01  1.85E+00  3.26E+01  3.22E+00 -2.05E-01  9.32E+01
 
 TH11
+        2.01E+09  2.37E+09  1.42E+09 -1.67E+09 -1.97E+09  1.85E+00  1.03E+01 -9.04E+08 -2.33E+09  6.35E+00  1.83E+09
 
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
 #CPUT: Total CPU Time in Seconds,       37.035
Stop Time:
Sat Sep 25 04:22:25 CDT 2021
