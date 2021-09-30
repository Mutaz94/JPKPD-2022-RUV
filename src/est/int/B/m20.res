Tue Sep 28 20:20:43 CDT 2021
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
$DATA ../../../../data/int/B/dat20.csv ignore=@
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
Current Date:       28 SEP 2021
Days until program expires : 201
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3320.92751098920        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0273E+02  4.9108E+01  7.5293E+01  7.3727E+01  1.2000E+02  3.6594E+01  1.2349E+01 -2.3973E+02 -4.0465E+01  7.4866E+00
            -7.7672E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3643.32710583693        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:      115
 NPARAMETR:  6.9234E-01  9.6041E-01  9.4992E-01  9.5725E-01  9.1103E-01  1.0319E+00  9.8262E-01  1.1852E+00  1.0268E+00  9.8631E-01
             1.4265E+00
 PARAMETER: -2.6768E-01  5.9604E-02  4.8623E-02  5.6307E-02  6.8154E-03  1.3144E-01  8.2466E-02  2.6990E-01  1.2646E-01  8.6212E-02
             4.5523E-01
 GRADIENT:  -6.9794E+02 -4.7639E+01 -1.6075E+01 -9.5471E+01 -3.6427E+01 -2.3047E+02  1.2862E+01  2.9326E+01 -5.1591E-02  1.5521E+01
             5.9094E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3670.24528297911        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  7.1806E-01  8.6299E-01  7.5973E-01  1.0412E+00  7.9247E-01  1.0196E+00  7.5679E-01  3.7895E-01  1.0260E+00  9.3628E-01
             1.4260E+00
 PARAMETER: -2.3121E-01 -4.7353E-02 -1.7479E-01  1.4033E-01 -1.3260E-01  1.1939E-01 -1.7866E-01 -8.7036E-01  1.2568E-01  3.4160E-02
             4.5487E-01
 GRADIENT:  -6.3944E+02 -2.6350E+01 -9.7350E+01 -1.9106E+01  2.2953E+01 -1.9155E+02 -3.7996E+00  2.6956E+00 -8.0105E+00 -1.9638E+00
             5.6168E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3832.27947280868        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      475
 NPARAMETR:  8.3497E-01  8.7838E-01  8.6454E-01  1.0157E+00  8.5616E-01  1.1170E+00  1.5299E+00  3.7638E-01  9.2955E-01  8.5625E-01
             1.1188E+00
 PARAMETER: -8.0359E-02 -2.9680E-02 -4.5559E-02  1.1557E-01 -5.5298E-02  2.1066E-01  5.2522E-01 -8.7716E-01  2.6943E-02 -5.5194E-02
             2.1227E-01
 GRADIENT:  -2.6646E+02 -2.8361E+01 -4.2728E+00 -2.7913E+01  5.0898E+00 -2.7600E+00  6.6844E+01 -5.0804E+00 -3.0634E+01 -7.0008E+00
             2.2257E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3852.16805581858        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      650
 NPARAMETR:  8.6357E-01  8.8991E-01  8.7209E-01  1.0166E+00  8.6145E-01  1.0970E+00  1.3980E+00  4.1738E-01  9.4168E-01  8.7193E-01
             1.0924E+00
 PARAMETER: -4.6676E-02 -1.6638E-02 -3.6868E-02  1.1649E-01 -4.9134E-02  1.9256E-01  4.3504E-01 -7.7375E-01  3.9910E-02 -3.7051E-02
             1.8839E-01
 GRADIENT:  -2.1417E+02 -2.6532E+01 -2.1793E+00 -2.1971E+01  4.0490E+00  5.7069E+00  5.3527E+01 -6.5516E+00 -2.6618E+01 -5.5735E+00
             1.8533E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3877.01641845003        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:      851            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7121E-01  9.1964E-01  8.7486E-01  1.0272E+00  8.6079E-01  1.0878E+00  9.9109E-01  4.1495E-01  1.0206E+00  8.9990E-01
             1.0940E+00
 PARAMETER:  7.0792E-02  1.6230E-02 -3.3691E-02  1.2688E-01 -4.9906E-02  1.8416E-01  9.1050E-02 -7.7959E-01  1.2037E-01 -5.4730E-03
             1.8981E-01
 GRADIENT:   3.9625E+02  6.2356E+01  1.0035E+01  1.1909E+02  3.8467E+01  1.5738E+02  1.2727E+01 -5.7613E+00  1.1160E+01 -5.8517E+00
             1.8556E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3878.97677259997        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1011
 NPARAMETR:  9.7222E-01  9.4231E-01  8.8222E-01  1.0204E+00  8.7471E-01  1.0101E+00  8.7636E-01  4.1491E-01  1.0470E+00  9.7514E-01
             1.0940E+00
 PARAMETER:  7.1825E-02  4.0576E-02 -2.5311E-02  1.2021E-01 -3.3867E-02  1.1000E-01 -3.1977E-02 -7.7968E-01  1.4596E-01  7.4823E-02
             1.8983E-01
 GRADIENT:   3.3843E+00  4.8511E+00 -2.8138E+00  6.8461E+00 -1.0645E+01 -1.0203E+00 -6.7668E-01 -6.3321E+00  5.4667E-01  1.1485E-01
             1.7833E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3883.91073173828        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1170
 NPARAMETR:  9.7022E-01  9.5253E-01  8.9457E-01  1.0107E+00  8.9124E-01  1.0117E+00  8.8186E-01  7.6616E-01  1.0447E+00  9.8191E-01
             1.0840E+00
 PARAMETER:  6.9765E-02  5.1364E-02 -1.1416E-02  1.1067E-01 -1.5140E-02  1.1167E-01 -2.5721E-02 -1.6637E-01  1.4377E-01  8.1744E-02
             1.8067E-01
 GRADIENT:   4.0197E+02  6.8081E+01 -1.7959E+01  9.6620E+01  4.3075E+01  8.0522E+01  5.3300E+00  1.8470E+00  2.0646E+01  1.0260E+01
             1.9691E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3884.15002806992        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  9.7022E-01  9.5253E-01  9.2278E-01  1.0107E+00  8.9124E-01  1.0045E+00  8.8155E-01  7.3978E-01  1.0447E+00  9.8157E-01
             1.0840E+00
 PARAMETER:  6.9765E-02  5.1364E-02  1.9636E-02  1.1067E-01 -1.5140E-02  1.0449E-01 -2.6069E-02 -2.0141E-01  1.4377E-01  8.1397E-02
             1.8067E-01
 GRADIENT:  -1.1685E+00 -5.0272E+00  7.7307E+00 -4.1326E+00 -2.3640E+01 -3.2980E+00  5.7065E-02 -1.9477E+00  3.2238E+00  5.1887E+00
             1.9123E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3884.34938219864        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1497
 NPARAMETR:  9.7125E-01  9.5549E-01  9.1726E-01  1.0127E+00  9.0280E-01  1.0139E+00  8.8030E-01  7.5825E-01  1.0444E+00  9.5571E-01
             1.0845E+00
 PARAMETER:  7.0833E-02  5.4469E-02  1.3639E-02  1.1259E-01 -2.2529E-03  1.1384E-01 -2.7498E-02 -1.7674E-01  1.4341E-01  5.4694E-02
             1.8112E-01
 GRADIENT:   1.1116E+00 -7.9985E+00 -7.5795E+00  4.1706E+00  8.7448E-01  4.6964E-01 -1.2688E+00 -1.0835E+00  2.8672E+00 -7.7992E-01
             1.9308E+02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -3884.43281318655        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1595
 NPARAMETR:  9.7092E-01  9.6171E-01  9.2559E-01  1.0095E+00  9.0454E-01  1.0130E+00  8.9327E-01  7.5819E-01  1.0443E+00  9.5960E-01
             1.0846E+00
 PARAMETER:  7.0406E-02  6.0888E-02  2.2695E-02  1.0967E-01 -1.3263E-03  1.1267E-01 -1.2353E-02 -1.7674E-01  1.4341E-01  5.8808E-02
             1.8112E-01
 GRADIENT:  -9.9191E-01 -5.9154E-01  1.5456E-01  1.6818E+00 -9.2657E+00 -5.1793E-01  2.7365E-01  1.6060E+05  9.9020E+04  1.3652E+04
            -1.5702E+05
 NUMSIGDIG:         2.0         2.1         2.8         1.8         1.0         1.6         1.3         2.3         2.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1595
 NO. OF SIG. DIGITS IN FINAL EST.:  1.0
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0814E-03 -2.6128E-02 -1.6380E-02  1.0030E-02 -1.8450E-02
 SE:             2.9926E-02  2.1467E-02  1.5768E-02  2.7883E-02  2.5071E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7117E-01  2.2355E-01  2.9892E-01  7.1907E-01  4.6179E-01

 ETASHRINKSD(%)  1.0000E-10  2.8083E+01  4.7174E+01  6.5889E+00  1.6010E+01
 ETASHRINKVR(%)  1.0000E-10  4.8279E+01  7.2094E+01  1.2744E+01  2.9457E+01
 EBVSHRINKSD(%)  2.9979E-01  2.7711E+01  4.9242E+01  6.2327E+00  1.5834E+01
 EBVSHRINKVR(%)  5.9869E-01  4.7743E+01  7.4237E+01  1.2077E+01  2.9160E+01
 RELATIVEINF(%)  9.9398E+01  2.4581E+01  1.5962E+01  6.3786E+01  2.7556E+01
 EPSSHRINKSD(%)  2.7824E+01
 EPSSHRINKVR(%)  4.7906E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3884.4328131865523     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2230.3434534181415     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.82
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3884.433       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.71E-01  9.62E-01  9.26E-01  1.01E+00  9.04E-01  1.01E+00  8.94E-01  7.58E-01  1.04E+00  9.60E-01  1.08E+00
 


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
+        1.14E+03
 
 TH 2
+       -2.24E+00  7.93E+02
 
 TH 3
+       -3.03E-01 -2.00E+01  5.61E+02
 
 TH 4
+        6.61E+07  2.58E+02  6.93E+07  8.75E+02
 
 TH 5
+       -1.38E+00 -4.31E+03 -3.73E+02 -1.76E+04  1.05E+03
 
 TH 6
+        2.82E+00 -1.35E+00  7.04E-01  5.62E+07 -9.29E-01  1.93E+02
 
 TH 7
+        4.52E-01  1.19E+01 -1.36E+01 -7.18E+07 -8.02E+00 -1.05E-01  8.89E+07
 
 TH 8
+        1.16E+03 -5.51E+07  2.12E+04 -4.78E+07  4.32E+01  1.89E+03 -4.44E+02  3.95E+07
 
 TH 9
+        1.04E+03 -4.93E+07  1.91E+04 -4.28E+07  5.25E+07  1.70E+03 -3.84E+02 -5.90E+04  3.16E+07
 
 TH10
+        4.46E-01  7.70E+07 -8.63E+00  6.68E+07 -2.90E+01  1.81E-01  2.73E+01  3.18E+00 -3.92E+00  1.12E+02
 
 TH11
+       -8.00E+02  1.74E+03 -1.46E+04  8.16E+03 -4.01E+07 -1.29E+03  3.23E+02  4.52E+04  1.46E+04  9.46E+00  1.85E+07
 
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
 #CPUT: Total CPU Time in Seconds,       56.016
Stop Time:
Tue Sep 28 20:21:44 CDT 2021
