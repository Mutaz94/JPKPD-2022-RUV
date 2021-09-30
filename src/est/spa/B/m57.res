Wed Sep 29 11:21:46 CDT 2021
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
$DATA ../../../../data/spa/B/dat57.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1716.59696031860        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6288E+02 -5.6660E+01 -5.2922E+01  3.8799E+00  6.4764E+01  6.9157E+01  1.1188E+01  1.2229E+01  3.1230E+01  2.5536E+01
            -3.3756E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1729.90566485284        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0252E+00  1.1150E+00  1.1479E+00  9.9655E-01  1.0415E+00  8.6962E-01  9.2246E-01  9.3826E-01  8.5902E-01  8.4280E-01
             1.1657E+00
 PARAMETER:  1.2494E-01  2.0889E-01  2.3791E-01  9.6547E-02  1.4061E-01 -3.9701E-02  1.9289E-02  3.6277E-02 -5.1964E-02 -7.1020E-02
             2.5331E-01
 GRADIENT:  -2.4244E+01  9.7075E+00  9.0171E+00  1.4321E+00  1.2319E-01 -2.0516E+01  1.5481E+00 -2.4071E+00 -7.0576E+00 -2.9158E+00
             1.9964E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1730.73226693635        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0358E+00  1.2460E+00  9.7364E-01  9.1733E-01  1.0079E+00  8.9081E-01  8.1781E-01  8.1248E-01  9.5203E-01  7.7359E-01
             1.1539E+00
 PARAMETER:  1.3518E-01  3.1993E-01  7.3287E-02  1.3712E-02  1.0782E-01 -1.5623E-02 -1.0112E-01 -1.0767E-01  5.0837E-02 -1.5672E-01
             2.4312E-01
 GRADIENT:   1.9399E+00  3.6561E+01  2.2089E+01  1.5316E+01 -2.9443E+01 -1.1189E+01 -1.3564E+00 -2.6376E+00 -3.8740E+00 -8.3688E+00
             1.4508E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1733.04582636459        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0376E+00  1.3701E+00  6.8365E-01  8.2275E-01  9.4097E-01  9.2158E-01  8.0422E-01  4.0192E-01  1.0149E+00  7.7237E-01
             1.1029E+00
 PARAMETER:  1.3695E-01  4.1487E-01 -2.8031E-01 -9.5101E-02  3.9158E-02  1.8333E-02 -1.1789E-01 -8.1150E-01  1.1481E-01 -1.5829E-01
             1.9791E-01
 GRADIENT:   2.2752E+00  2.5583E+01  3.7827E+00  1.7994E+01 -1.9752E+01  9.6578E-01  7.4506E-01  6.4053E-01  1.3008E+00  3.5470E+00
             3.5471E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1734.11667145443        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0351E+00  1.7052E+00  5.3338E-01  5.9975E-01  1.0675E+00  9.2020E-01  6.6030E-01  1.3601E-01  1.2800E+00  8.2136E-01
             1.0871E+00
 PARAMETER:  1.3453E-01  6.3370E-01 -5.2852E-01 -4.1124E-01  1.6530E-01  1.6832E-02 -3.1506E-01 -1.8950E+00  3.4686E-01 -9.6798E-02
             1.8350E-01
 GRADIENT:  -4.4413E+00  1.2139E+01 -1.3555E+00  1.1557E+01  1.6049E+00 -1.6859E-01 -3.2099E+00  4.9613E-02 -1.3904E+00 -1.3516E+00
            -4.6422E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1734.22693096969        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0361E+00  1.8288E+00  5.0055E-01  5.1859E-01  1.1384E+00  9.2123E-01  6.2600E-01  7.7786E-02  1.4411E+00  8.7472E-01
             1.0974E+00
 PARAMETER:  1.3544E-01  7.0365E-01 -5.9204E-01 -5.5665E-01  2.2965E-01  1.7951E-02 -3.6840E-01 -2.4538E+00  4.6538E-01 -3.3852E-02
             1.9297E-01
 GRADIENT:  -1.9276E+00  5.5478E+00 -1.5584E+00  8.9884E+00  6.4942E+00  3.1559E-01 -3.0658E+00  1.2405E-02 -3.5874E-01 -1.0568E+00
            -1.3916E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1734.33111156474        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1083             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0378E+00  1.8447E+00  4.9024E-01  4.9905E-01  1.1448E+00  9.2042E-01  6.3456E-01  1.3613E-02  1.4682E+00  8.7909E-01
             1.0998E+00
 PARAMETER:  1.3709E-01  7.1232E-01 -6.1286E-01 -5.9506E-01  2.3526E-01  1.7077E-02 -3.5482E-01 -4.1967E+00  4.8406E-01 -2.8865E-02
             1.9516E-01
 GRADIENT:   5.0472E+02  6.6816E+02  2.8368E+00  8.5001E+01  1.4395E+01  3.1036E+01  1.2658E+01  7.2600E-04  1.4677E+01  3.9326E-01
             1.5843E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1734.33385679681        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1269             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0377E+00  1.8457E+00  4.8920E-01  4.9937E-01  1.1430E+00  9.2040E-01  6.3467E-01  1.0000E-02  1.4670E+00  8.7761E-01
             1.0996E+00
 PARAMETER:  1.3705E-01  7.1286E-01 -6.1498E-01 -5.9440E-01  2.3362E-01  1.7054E-02 -3.5465E-01 -5.2233E+00  4.8322E-01 -3.0556E-02
             1.9494E-01
 GRADIENT:   5.0458E+02  6.7230E+02  3.3068E+00  8.5514E+01  1.2467E+01  3.1038E+01  1.2672E+01  0.0000E+00  1.4582E+01  3.9925E-01
             1.4367E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1734.33385679681        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:     1330
 NPARAMETR:  1.0377E+00  1.8437E+00  4.8620E-01  4.9986E-01  1.1428E+00  9.2041E-01  6.3551E-01  1.0000E-02  1.4677E+00  8.7814E-01
             1.1003E+00
 PARAMETER:  1.3705E-01  7.1286E-01 -6.1498E-01 -5.9440E-01  2.3362E-01  1.7054E-02 -3.5465E-01 -5.2233E+00  4.8322E-01 -3.0556E-02
             1.9494E-01
 GRADIENT:  -3.2036E-03  1.6123E+00  4.9803E-01 -2.4652E-01  1.0634E-01 -1.4984E-03 -8.1904E-02  0.0000E+00 -3.3069E-02 -3.6347E-02
            -1.2156E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1330
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8153E-04 -3.3026E-02 -2.2954E-04  2.6819E-02 -3.6035E-02
 SE:             2.9800E-02  2.2813E-02  9.2462E-05  2.2822E-02  2.1876E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9246E-01  1.4771E-01  1.3047E-02  2.3994E-01  9.9518E-02

 ETASHRINKSD(%)  1.6545E-01  2.3572E+01  9.9690E+01  2.3542E+01  2.6711E+01
 ETASHRINKVR(%)  3.3063E-01  4.1588E+01  9.9999E+01  4.1542E+01  4.6288E+01
 EBVSHRINKSD(%)  5.9840E-01  2.2899E+01  9.9721E+01  2.5342E+01  2.5780E+01
 EBVSHRINKVR(%)  1.1932E+00  4.0554E+01  9.9999E+01  4.4262E+01  4.4914E+01
 RELATIVEINF(%)  9.8749E+01  3.6755E+00  8.3019E-05  3.6450E+00  1.1300E+01
 EPSSHRINKSD(%)  4.2742E+01
 EPSSHRINKVR(%)  6.7216E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1734.3338567968131     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -999.18303023307487     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.58
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1734.334       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.85E+00  4.89E-01  4.99E-01  1.14E+00  9.20E-01  6.35E-01  1.00E-02  1.47E+00  8.78E-01  1.10E+00
 


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
+        1.20E+03
 
 TH 2
+       -1.04E+01  5.04E+02
 
 TH 3
+        9.59E+00  1.85E+02  3.89E+02
 
 TH 4
+       -2.46E+01  4.23E+02 -3.10E+02  1.17E+03
 
 TH 5
+       -5.57E+00 -2.43E+02 -3.65E+02  2.74E+02  6.09E+02
 
 TH 6
+        5.62E-02 -1.49E+00  2.36E+00 -5.93E+00 -1.23E+00  2.29E+02
 
 TH 7
+        1.38E+00 -2.30E+00 -1.91E+00 -1.50E+01 -2.27E+01  9.28E-01  1.77E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.29E+00 -1.92E+01 -3.11E+01  6.16E+01  5.98E-01 -5.09E-01  2.41E+01  0.00E+00  3.79E+01
 
 TH10
+       -1.58E-01 -1.59E+01 -3.87E+01 -3.28E+00 -6.48E+01  3.75E-02  2.46E+01  0.00E+00  4.05E+00  9.02E+01
 
 TH11
+       -8.50E+00 -2.25E+01 -2.51E+01  1.40E+00 -4.73E+00  3.15E+00  1.43E+01  0.00E+00  5.52E+00  2.16E+01  1.81E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.641
Stop Time:
Wed Sep 29 11:22:11 CDT 2021
