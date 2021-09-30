Wed Sep 29 12:58:04 CDT 2021
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
$DATA ../../../../data/spa/A2/dat65.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -690.065791514912        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0918E+02  7.1672E+01  2.4749E+01  1.2283E+02  2.8672E+02  7.7855E+01 -5.6003E+01 -3.9833E+01 -5.5682E+01 -1.6520E+02
            -1.6630E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1308.67950913850        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0873E+00  8.6959E-01  8.1366E-01  1.0686E+00  7.6608E-01  8.9632E-01  1.1393E+00  1.0489E+00  1.2689E+00  1.1828E+00
             1.9791E+00
 PARAMETER:  1.8371E-01 -3.9730E-02 -1.0621E-01  1.6633E-01 -1.6647E-01 -9.4587E-03  2.3045E-01  1.4776E-01  3.3814E-01  2.6793E-01
             7.8264E-01
 GRADIENT:   3.5447E+02  1.1212E+01 -5.4036E+00  5.2517E+01  5.3888E+01  1.6315E+01 -3.2926E+00  3.9407E+00  3.5198E+01  2.1081E+00
            -2.5668E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1317.23892264619        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      204
 NPARAMETR:  1.0900E+00  7.4313E-01  5.1379E-01  1.1273E+00  5.2564E-01  8.1087E-01  1.6244E+00  6.6361E-01  1.0820E+00  9.1411E-01
             1.9262E+00
 PARAMETER:  1.8618E-01 -1.9689E-01 -5.6595E-01  2.1986E-01 -5.4315E-01 -1.0965E-01  5.8511E-01 -3.1006E-01  1.7877E-01  1.0193E-02
             7.5554E-01
 GRADIENT:   1.6572E+02  4.6568E+01  1.7996E+00  5.0919E+01  2.0011E+01 -3.0801E+01  1.4033E+01  2.7679E+00  4.8482E+00  2.5081E+00
            -2.6554E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1370.66057585620        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  1.0627E+00  6.5380E-01  3.2695E-01  1.0779E+00  3.8529E-01  8.6267E-01  9.6897E-01  3.8631E-01  1.0837E+00  6.3852E-01
             2.3940E+00
 PARAMETER:  1.6080E-01 -3.2496E-01 -1.0180E+00  1.7503E-01 -8.5377E-01 -4.7720E-02  6.8476E-02 -8.5111E-01  1.8036E-01 -3.4860E-01
             9.7296E-01
 GRADIENT:   4.9051E+01  3.8240E+01  1.9381E+01  4.2616E+01 -1.6946E+01 -1.3251E+00 -1.1580E+01 -1.1295E+00 -1.4054E+00 -7.3428E+00
            -8.5414E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1382.39233062836        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      556
 NPARAMETR:  1.0528E+00  4.3474E-01  3.0547E-01  1.1274E+00  3.2232E-01  8.5060E-01  1.3109E+00  1.2820E-01  1.0395E+00  6.6089E-01
             2.7114E+00
 PARAMETER:  1.5148E-01 -7.3302E-01 -1.0859E+00  2.1991E-01 -1.0322E+00 -6.1810E-02  3.7069E-01 -1.9542E+00  1.3878E-01 -3.1416E-01
             1.0975E+00
 GRADIENT:   5.0889E+00  1.9312E+00 -6.1502E+00  7.4179E+00  1.0894E+01 -3.4179E+00 -1.0126E+00 -9.3579E-02  2.2619E+00 -1.4498E+00
             3.2949E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1384.13994600025        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  1.0459E+00  2.2979E-01  3.6267E-01  1.2430E+00  3.1349E-01  8.5039E-01  1.8849E+00  5.5865E-02  9.9386E-01  7.6725E-01
             2.6957E+00
 PARAMETER:  1.4488E-01 -1.3706E+00 -9.1427E-01  3.1749E-01 -1.0600E+00 -6.2061E-02  7.3387E-01 -2.7848E+00  9.3839E-02 -1.6494E-01
             1.0917E+00
 GRADIENT:  -4.6966E-01  5.2253E+00  2.4417E+01  8.4463E-01 -3.6525E+01  1.1147E+00 -4.3603E-01 -2.0104E-03  2.8730E+00 -8.2360E-01
            -3.6598E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1385.58762201413        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.0369E+00  8.1908E-02  4.0042E-01  1.3324E+00  3.2591E-01  8.3789E-01  4.0311E+00  1.0000E-02  9.4630E-01  7.8464E-01
             2.7365E+00
 PARAMETER:  1.3619E-01 -2.4022E+00 -8.1523E-01  3.8701E-01 -1.0211E+00 -7.6868E-02  1.4940E+00 -4.5602E+00  4.4804E-02 -1.4253E-01
             1.1067E+00
 GRADIENT:   1.3622E+00  2.0477E-01 -4.1311E+00  7.8339E+00  5.0653E+00 -9.4221E-01 -2.5034E-01  0.0000E+00  1.1799E+00 -5.2885E-01
             7.5051E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1385.68822982366        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  1.0345E+00  5.6722E-02  3.9487E-01  1.3323E+00  3.1892E-01  8.3985E-01  5.0341E+00  1.0000E-02  9.4030E-01  7.9246E-01
             2.7275E+00
 PARAMETER:  1.3392E-01 -2.7696E+00 -8.2919E-01  3.8693E-01 -1.0428E+00 -7.4528E-02  1.7162E+00 -5.2381E+00  3.8445E-02 -1.3262E-01
             1.1034E+00
 GRADIENT:  -1.2438E+00  3.8214E-01  1.8905E+00  1.3834E+00 -3.7450E+00  8.9547E-02  3.0150E-01  0.0000E+00 -1.1131E-01  1.1685E-01
            -7.5892E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1385.73193664724        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1259
 NPARAMETR:  1.0322E+00  1.9512E-02  4.0648E-01  1.3557E+00  3.2209E-01  8.3813E-01  8.5913E+00  1.0000E-02  9.3231E-01  7.9883E-01
             2.7315E+00
 PARAMETER:  1.3165E-01 -3.8367E+00 -8.0022E-01  4.0435E-01 -1.0329E+00 -7.6576E-02  2.2508E+00 -6.9716E+00  2.9913E-02 -1.2461E-01
             1.1049E+00
 GRADIENT:  -7.3935E-01  8.6675E-02  7.3910E-01  2.5899E+00 -1.6837E+00 -7.0774E-02  3.7275E-02  0.0000E+00  4.4020E-03  5.5823E-03
            -1.7212E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1385.74196300740        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1429
 NPARAMETR:  1.0322E+00  1.0000E-02  4.0698E-01  1.3560E+00  3.2222E-01  8.3817E-01  1.1699E+01  1.0000E-02  9.3109E-01  7.9976E-01
             2.7321E+00
 PARAMETER:  1.3112E-01 -4.5224E+00 -7.9757E-01  4.0703E-01 -1.0339E+00 -7.6932E-02  2.5853E+00 -8.1111E+00  2.7839E-02 -1.2267E-01
             1.1053E+00
 GRADIENT:  -8.0907E-01  7.2348E-03  7.4698E-01  2.6945E+00 -1.3979E+00 -6.7165E-02  3.2653E+00  0.0000E+00 -1.0591E-01  5.0135E-02
             7.0372E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1429
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.7878E-04  1.5099E-03  3.6734E-05 -1.0741E-02 -3.9789E-03
 SE:             2.8717E-02  1.8443E-03  2.1015E-04  2.6872E-02  2.2455E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8948E-01  4.1295E-01  8.6124E-01  6.8936E-01  8.5936E-01

 ETASHRINKSD(%)  3.7960E+00  9.3822E+01  9.9296E+01  9.9763E+00  2.4773E+01
 ETASHRINKVR(%)  7.4479E+00  9.9618E+01  9.9995E+01  1.8957E+01  4.3409E+01
 EBVSHRINKSD(%)  3.6786E+00  9.4612E+01  9.9307E+01  8.9802E+00  2.3919E+01
 EBVSHRINKVR(%)  7.2219E+00  9.9710E+01  9.9995E+01  1.7154E+01  4.2116E+01
 RELATIVEINF(%)  7.0237E+01  2.0781E-02  2.1385E-04  9.1743E+00  2.2418E+00
 EPSSHRINKSD(%)  3.2919E+01
 EPSSHRINKVR(%)  5.5001E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1385.7419630073962     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -650.59113644365800     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1385.742       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.00E-02  4.08E-01  1.36E+00  3.22E-01  8.38E-01  1.20E+01  1.00E-02  9.30E-01  8.00E-01  2.73E+00
 


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
+        1.38E+03
 
 TH 2
+       -3.18E+01  9.97E+04
 
 TH 3
+       -5.71E+00  1.17E+02  3.19E+03
 
 TH 4
+       -5.43E+01  5.24E+01 -3.34E+02  5.79E+02
 
 TH 5
+        1.53E+02 -2.71E+02 -5.09E+03 -2.29E+02  9.76E+03
 
 TH 6
+       -2.73E+00 -2.07E+00  1.90E+01 -1.28E+01 -9.89E-01  2.42E+02
 
 TH 7
+       -6.76E-01 -1.45E+01  1.22E+01 -1.02E+01  7.00E+00 -9.86E-02  1.87E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.12E+01 -1.33E+01  5.34E+01 -1.25E+01  1.15E+01  5.18E+00 -3.68E-01  0.00E+00  1.61E+02
 
 TH10
+       -9.49E+00 -1.33E+01 -5.09E+01  6.66E+00  2.03E+01 -3.58E-01  5.04E+00  0.00E+00  5.19E-01  1.01E+02
 
 TH11
+       -1.82E+01 -9.37E+00 -1.24E+01 -8.17E+00  7.61E+00  4.09E+00 -1.69E+00  0.00E+00  6.22E+00  2.04E+01  3.84E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.653
Stop Time:
Wed Sep 29 12:58:31 CDT 2021
