Wed Sep 29 03:34:12 CDT 2021
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
$DATA ../../../../data/int/SL2/dat80.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      996
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

 TOT. NO. OF OBS RECS:      896
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2125.83239426925        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1836E+02  1.3657E+02  1.2503E+02  1.0614E+02  9.2329E+01  8.4432E+01 -6.2328E+01 -1.2836E+02 -2.3840E+01  1.7368E+00
            -3.2771E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3066.27497111962        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0903E+00  1.0758E+00  9.4641E-01  1.0310E+00  9.9588E-01  8.2838E-01  9.9590E-01  1.0044E+00  8.9634E-01  9.1004E-01
             2.2564E+00
 PARAMETER:  1.8643E-01  1.7303E-01  4.4925E-02  1.3050E-01  9.5872E-02 -8.8283E-02  9.5890E-02  1.0434E-01 -9.4329E-03  5.7352E-03
             9.1377E-01
 GRADIENT:   3.4190E+02  8.5727E+01 -9.1952E+00  1.0655E+02 -1.5683E+01 -2.7933E+01  3.6540E+00  6.3058E+00  2.2976E+00  7.9394E-01
             1.6892E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3069.33007799722        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0754E+00  1.0171E+00  8.5107E-01  1.0492E+00  9.3034E-01  8.4645E-01  1.0366E+00  7.4175E-01  9.2142E-01  7.7836E-01
             2.2320E+00
 PARAMETER:  1.7274E-01  1.1693E-01 -6.1264E-02  1.4807E-01  2.7797E-02 -6.6709E-02  1.3591E-01 -1.9875E-01  1.8166E-02 -1.5057E-01
             9.0290E-01
 GRADIENT:   2.8473E+02  6.7013E+01 -2.2058E+01  9.9507E+01  4.0410E+00 -1.5042E+01 -2.4726E+00  1.7770E+00  8.0871E+00 -1.1130E+01
             1.4428E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3079.42803204723        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      248
 NPARAMETR:  1.0068E+00  1.0082E+00  9.6165E-01  1.0057E+00  9.8197E-01  8.5122E-01  1.0386E+00  8.8065E-01  9.0979E-01  8.4497E-01
             2.1017E+00
 PARAMETER:  1.0676E-01  1.0813E-01  6.0891E-02  1.0571E-01  8.1806E-02 -6.1087E-02  1.3785E-01 -2.7090E-02  5.4628E-03 -6.8459E-02
             8.4274E-01
 GRADIENT:  -6.8157E+01 -7.2120E+00  2.9008E+00 -1.2351E+01 -8.6937E+00 -1.7653E+01 -6.5779E-01 -8.2039E-04  8.3455E-01 -6.7137E+00
             1.8902E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3083.34928149780        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      427
 NPARAMETR:  1.0345E+00  1.3485E+00  1.2567E+00  8.3594E-01  1.3432E+00  8.9025E-01  8.1184E-01  1.5462E+00  9.8972E-01  1.0917E+00
             2.0836E+00
 PARAMETER:  1.3388E-01  3.9901E-01  3.2852E-01 -7.9197E-02  3.9504E-01 -1.6249E-02 -1.0845E-01  5.3580E-01  8.9665E-02  1.8775E-01
             8.3408E-01
 GRADIENT:   1.1903E+01  1.1752E+01 -4.4621E+00  1.6614E+01  7.7934E+00  1.4530E+00 -1.0264E-01 -1.5867E-01 -1.0317E+00 -3.4601E+00
            -4.3175E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3083.75481137245        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  1.0328E+00  1.3952E+00  1.5033E+00  7.9872E-01  1.4368E+00  8.8882E-01  7.8077E-01  2.0674E+00  1.0005E+00  1.1527E+00
             2.0813E+00
 PARAMETER:  1.3229E-01  4.3302E-01  5.0767E-01 -1.2475E-01  4.6239E-01 -1.7856E-02 -1.4748E-01  8.2629E-01  1.0055E-01  2.4211E-01
             8.3299E-01
 GRADIENT:   8.4727E+00 -5.2740E+00 -1.1159E+00 -7.2111E+00 -2.5217E+00  1.1199E+00 -1.0639E+00  8.5025E-01  1.3632E-01  2.7274E+00
             1.7703E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3085.01779886023        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      787             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0303E+00  1.2964E+00  1.9864E+00  8.7414E-01  1.4543E+00  8.8577E-01  8.5087E-01  2.6509E+00  8.9770E-01  1.0762E+00
             2.0668E+00
 PARAMETER:  1.2980E-01  3.5956E-01  7.8634E-01 -3.4516E-02  4.7455E-01 -2.1294E-02 -6.1501E-02  1.0749E+00 -7.9236E-03  1.7346E-01
             8.2602E-01
 GRADIENT:   1.4739E+02  8.6858E+01  4.0990E+00  1.4112E+01  4.1037E+01  1.0431E+01  1.0199E+00  4.9734E+00  1.7269E+00  1.6258E+00
             1.6742E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3085.04911284274        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      966
 NPARAMETR:  1.0295E+00  1.2880E+00  1.9864E+00  8.8103E-01  1.4544E+00  8.8611E-01  8.6360E-01  2.6508E+00  8.8728E-01  1.0722E+00
             2.0644E+00
 PARAMETER:  1.2906E-01  3.5306E-01  7.8630E-01 -2.6663E-02  4.7457E-01 -2.0919E-02 -4.6648E-02  1.0748E+00 -1.9591E-02  1.6974E-01
             8.2486E-01
 GRADIENT:   2.8533E-01  4.3567E-01 -3.8029E+00 -9.1409E-01  5.4579E+00  3.0891E-02 -1.1291E-02 -6.9800E-02  3.4658E-03  3.7155E-02
            -2.4263E-01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -3085.04911284274        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:      990
 NPARAMETR:  1.0295E+00  1.2880E+00  1.9864E+00  8.8103E-01  1.4544E+00  8.8611E-01  8.6360E-01  2.6508E+00  8.8728E-01  1.0722E+00
             2.0644E+00
 PARAMETER:  1.2906E-01  3.5306E-01  7.8630E-01 -2.6663E-02  4.7457E-01 -2.0919E-02 -4.6648E-02  1.0748E+00 -1.9591E-02  1.6974E-01
             8.2486E-01
 GRADIENT:   2.6185E+04 -9.5814E+03  4.3036E+03 -5.2462E-01 -7.1175E+03  2.3883E-04  3.3799E+04  3.0834E+03  1.3250E-03 -1.9911E+04
            -4.1057E+03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      990
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6733E-03 -2.3554E-02 -2.7996E-02  1.6492E-02 -3.8886E-02
 SE:             2.9617E-02  2.1011E-02  1.9891E-02  2.2930E-02  2.2139E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5494E-01  2.6227E-01  1.5930E-01  4.7200E-01  7.9005E-02

 ETASHRINKSD(%)  7.7826E-01  2.9612E+01  3.3362E+01  2.3183E+01  2.5833E+01
 ETASHRINKVR(%)  1.5505E+00  5.0455E+01  5.5594E+01  4.0991E+01  4.4993E+01
 EBVSHRINKSD(%)  1.1921E+00  2.9814E+01  3.7211E+01  2.5681E+01  2.2923E+01
 EBVSHRINKVR(%)  2.3700E+00  5.0739E+01  6.0575E+01  4.4767E+01  4.0591E+01
 RELATIVEINF(%)  9.7587E+01  8.5764E+00  2.1520E+01  1.0383E+01  2.4293E+01
 EPSSHRINKSD(%)  1.8428E+01
 EPSSHRINKVR(%)  3.3460E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          896
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1646.7378515027735     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3085.0491128427443     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1438.3112613399708     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3085.049       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.29E+00  1.99E+00  8.81E-01  1.45E+00  8.86E-01  8.64E-01  2.65E+00  8.87E-01  1.07E+00  2.06E+00
 


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
+        4.79E+06
 
 TH 2
+       -7.86E+01  4.10E+05
 
 TH 3
+        2.00E+01 -1.19E+05  3.48E+04
 
 TH 4
+       -7.22E+06  1.90E+03 -4.37E+02  1.09E+07
 
 TH 5
+       -4.79E+01  2.69E+05 -7.78E+01  1.01E+03  1.78E+05
 
 TH 6
+       -7.18E+06 -7.18E+01  2.02E+01  1.08E+07 -4.62E+01  2.42E+02
 
 TH 7
+        7.36E+06 -5.11E+02  1.53E+02 -1.11E+07 -3.50E+02  3.61E+02  1.13E+07
 
 TH 8
+        1.06E+01 -1.36E+02  7.71E+01 -2.22E+02 -2.86E+01  1.09E+01  8.49E+01  1.00E+04
 
 TH 9
+        7.17E+06 -1.10E+02  2.97E+01 -1.08E+07 -5.73E+01 -2.95E-01  5.63E+02  2.01E+01  8.26E+01
 
 TH10
+       -1.65E+02  1.07E+02 -3.56E+01  3.49E+03  4.74E+01 -1.71E+02 -1.31E+03 -1.62E+01 -2.50E+02  2.55E+06
 
 TH11
+       -3.25E+01  4.27E+02 -1.30E+02  3.67E+02  7.20E+04 -1.54E+01 -1.38E+02 -6.45E+01 -1.81E+01  3.97E+01  2.95E+04
 
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
 #CPUT: Total CPU Time in Seconds,       41.303
Stop Time:
Wed Sep 29 03:34:55 CDT 2021
