Wed Sep 29 14:13:34 CDT 2021
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
$DATA ../../../../data/spa/S1/dat37.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1636.67121017797        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3605E+02  1.6073E+01  8.7555E+00  5.2161E+01  1.8615E+01  6.1223E+01 -2.5927E+01 -6.0758E+00 -1.0478E+01 -1.4787E+01
            -1.4307E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1642.66691095865        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      196
 NPARAMETR:  9.9747E-01  1.0897E+00  9.2683E-01  9.2796E-01  1.0340E+00  9.7835E-01  1.4448E+00  1.0442E+00  1.0715E+00  1.1003E+00
             1.0309E+00
 PARAMETER:  9.7462E-02  1.8594E-01  2.4011E-02  2.5234E-02  1.3345E-01  7.8116E-02  4.6794E-01  1.4327E-01  1.6908E-01  1.9557E-01
             1.3042E-01
 GRADIENT:  -2.4121E+01 -2.1212E+01 -9.5700E+00 -2.4797E+01  7.0652E+00 -3.6816E+00  1.2186E+01  4.8185E+00  1.8802E+01  8.2643E+00
             5.3531E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1645.20414975560        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.9525E-01  8.6196E-01  1.0104E+00  1.0867E+00  9.6151E-01  9.5929E-01  1.7070E+00  1.0215E+00  8.8182E-01  1.1158E+00
             1.0337E+00
 PARAMETER:  9.5237E-02 -4.8546E-02  1.1036E-01  1.8317E-01  6.0749E-02  5.8439E-02  6.3475E-01  1.2130E-01 -2.5770E-02  2.0959E-01
             1.3310E-01
 GRADIENT:  -2.6369E+01 -4.5029E+00 -1.0315E+01 -1.3418E+01 -1.8543E+00 -1.0587E+01  4.7704E+00  5.0461E+00  9.0865E+00  1.4021E+01
             6.3399E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1647.34765236694        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      552
 NPARAMETR:  1.0093E+00  7.3618E-01  9.7337E-01  1.1686E+00  8.7770E-01  9.8393E-01  1.9146E+00  8.8589E-01  8.0684E-01  9.5308E-01
             1.0257E+00
 PARAMETER:  1.0924E-01 -2.0629E-01  7.3009E-02  2.5582E-01 -3.0455E-02  8.3802E-02  7.4949E-01 -2.1166E-02 -1.1463E-01  5.1942E-02
             1.2534E-01
 GRADIENT:   8.8083E+00  6.1192E+00  1.4695E+00  2.9581E+00 -3.6773E+00  2.2929E-01  1.6274E+00  1.0507E-01  1.4176E+00  3.4406E-01
            -5.9868E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1647.53094182652        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      735             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0022E+00  6.6351E-01  1.0184E+00  1.2136E+00  8.7358E-01  9.8252E-01  2.0273E+00  9.5490E-01  7.9523E-01  9.4879E-01
             1.0304E+00
 PARAMETER:  1.0222E-01 -3.1021E-01  1.1827E-01  2.9359E-01 -3.5154E-02  8.2362E-02  8.0672E-01  5.3852E-02 -1.2913E-01  4.7433E-02
             1.2999E-01
 GRADIENT:   4.1892E+02  5.2097E+01  2.6023E+00  3.5252E+02  5.0489E+00  4.8820E+01  4.5937E+01  1.5769E+00  1.3785E+01  4.4223E-01
             3.0351E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1647.57349401255        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      891
 NPARAMETR:  1.0022E+00  6.6351E-01  1.0184E+00  1.2136E+00  8.7358E-01  9.8216E-01  2.0273E+00  9.3570E-01  7.8196E-01  9.4879E-01
             1.0266E+00
 PARAMETER:  1.0222E-01 -3.1021E-01  1.1827E-01  2.9359E-01 -3.5154E-02  8.1999E-02  8.0672E-01  3.3536E-02 -1.4595E-01  4.7434E-02
             1.2623E-01
 GRADIENT:  -4.9545E+00  4.0852E+00  1.7364E+00  5.8055E+00 -3.1616E+00  1.2807E-02 -1.6545E+00  1.2562E-04  7.0944E-02 -1.0948E+00
            -3.6053E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1647.64798278709        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  1.0023E+00  6.5183E-01  1.0178E+00  1.2141E+00  8.7460E-01  9.8040E-01  2.0628E+00  9.2246E-01  7.7461E-01  9.6184E-01
             1.0256E+00
 PARAMETER:  1.0233E-01 -3.2798E-01  1.1768E-01  2.9403E-01 -3.3990E-02  8.0201E-02  8.2406E-01  1.9294E-02 -1.5539E-01  6.1096E-02
             1.2527E-01
 GRADIENT:  -4.3168E+00  7.4379E-01 -4.6967E-01 -3.3630E+00  7.4719E-01 -5.9353E-01 -6.7025E-01 -1.6048E-01 -4.1990E-01  5.7776E-02
            -1.6322E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1647.65382728711        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  9.9981E-01  6.3738E-01  1.0217E+00  1.2204E+00  8.7165E-01  9.7906E-01  2.0875E+00  9.2486E-01  7.7162E-01  9.6159E-01
             1.0257E+00
 PARAMETER:  9.9811E-02 -3.5039E-01  1.2151E-01  2.9920E-01 -3.7372E-02  7.8842E-02  8.3598E-01  2.1890E-02 -1.5926E-01  6.0837E-02
             1.2536E-01
 GRADIENT:  -9.7415E+00 -6.1380E-01 -1.3177E-01 -6.8186E+00  2.0327E-01 -1.0741E+00 -1.1490E+00 -1.6924E-01 -3.6844E-01  1.9024E-01
            -6.1215E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1647.70205608297        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1408
 NPARAMETR:  1.0050E+00  6.3328E-01  1.0232E+00  1.2259E+00  8.7088E-01  9.8142E-01  2.1061E+00  9.2908E-01  7.7243E-01  9.6020E-01
             1.0260E+00
 PARAMETER:  1.0498E-01 -3.5684E-01  1.2292E-01  3.0368E-01 -3.8248E-02  8.1247E-02  8.4483E-01  2.6441E-02 -1.5821E-01  5.9385E-02
             1.2568E-01
 GRADIENT:   2.4329E+00  7.9288E-01 -1.2521E+00 -1.6709E+00  1.2974E+00 -6.1843E-02 -5.0209E-01  1.0602E-02  1.5345E-01  4.7665E-02
             7.5688E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1647.73424658722        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1585
 NPARAMETR:  9.9835E-01  5.8655E-01  1.0430E+00  1.2562E+00  8.6350E-01  9.8000E-01  2.2118E+00  9.4366E-01  7.5955E-01  9.5989E-01
             1.0262E+00
 PARAMETER:  9.8348E-02 -4.3350E-01  1.4207E-01  3.2806E-01 -4.6765E-02  7.9793E-02  8.9379E-01  4.2007E-02 -1.7503E-01  5.9067E-02
             1.2583E-01
 GRADIENT:  -1.1286E+01  1.4200E+00 -7.8054E-01  3.0123E+00  5.3771E-03 -3.8984E-01 -1.3675E+00 -3.3074E-01 -4.5631E-01 -1.7220E-01
            -2.3772E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1647.73580031094        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1650
 NPARAMETR:  9.9803E-01  5.8398E-01  1.0435E+00  1.2565E+00  8.6344E-01  9.8098E-01  2.2125E+00  9.4439E-01  7.5939E-01  9.6009E-01
             1.0262E+00
 PARAMETER:  9.8248E-02 -4.3694E-01  1.4289E-01  3.2902E-01 -4.7050E-02  7.9801E-02  8.9609E-01  4.3007E-02 -1.7562E-01  5.9052E-02
             1.2584E-01
 GRADIENT:   2.4609E+05  2.8162E+04  1.7227E+05  3.7394E+04 -2.4610E+05 -4.4716E-01  2.7464E+04  2.4599E+05 -7.0072E+04 -1.2305E+05
            -2.2000E-02
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         1.6         2.3         2.3         2.3         2.3
                    3.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1650
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.7047E-03  2.2226E-02 -3.5388E-02 -2.9032E-02 -2.0318E-02
 SE:             2.9885E-02  2.1107E-02  1.5306E-02  2.2619E-02  2.0688E-02
 N:                     100         100         100         100         100

 P VAL.:         8.2248E-01  2.9233E-01  2.0779E-02  1.9932E-01  3.2603E-01

 ETASHRINKSD(%)  1.0000E-10  2.9290E+01  4.8721E+01  2.4223E+01  3.0694E+01
 ETASHRINKVR(%)  1.0000E-10  5.0001E+01  7.3705E+01  4.2578E+01  5.1966E+01
 EBVSHRINKSD(%)  4.9630E-01  3.1338E+01  5.2353E+01  2.3176E+01  2.7293E+01
 EBVSHRINKVR(%)  9.9014E-01  5.2856E+01  7.7298E+01  4.0981E+01  4.7137E+01
 RELATIVEINF(%)  9.8288E+01  6.8349E+00  4.0817E+00  8.9746E+00  9.4177E+00
 EPSSHRINKSD(%)  4.4986E+01
 EPSSHRINKVR(%)  6.9735E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1647.7358003109362     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -912.58497374719798     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.04
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1647.736       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.98E-01  5.85E-01  1.04E+00  1.26E+00  8.63E-01  9.80E-01  2.22E+00  9.45E-01  7.59E-01  9.60E-01  1.03E+00
 


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
+        6.17E+07
 
 TH 2
+        3.63E+02  9.43E+06
 
 TH 3
+        6.54E+02  7.12E+03  2.77E+07
 
 TH 4
+        2.29E+02 -2.52E+02  4.21E+03  3.59E+06
 
 TH 5
+       -1.11E+03  1.95E+03 -2.12E+04  1.47E+03  8.26E+07
 
 TH 6
+        1.19E+03  4.57E+02  7.90E+02  2.83E+02 -1.36E+03  2.04E+02
 
 TH 7
+        3.68E+02  1.75E+02  2.45E+02  7.59E+01 -4.22E+02  5.93E+01  1.56E+05
 
 TH 8
+       -6.11E+04 -2.39E+04 -4.10E+04 -1.47E+04  7.07E+04  1.24E+03 -3.07E+03  6.89E+07
 
 TH 9
+       -5.02E+03 -1.98E+03 -3.37E+03 -1.23E+03  5.82E+03 -8.80E+02 -2.43E+02 -4.88E+07  3.46E+07
 
 TH10
+       -3.45E+02 -1.33E+02 -2.45E+02 -1.01E+02  3.17E+02 -1.22E+03 -1.71E+01  6.36E+04  5.22E+03  6.68E+07
 
 TH11
+       -1.77E+05 -6.93E+04 -1.19E+05 -4.28E+04  2.05E+05  2.14E+00 -8.90E+03 -1.87E+05  1.33E+05  1.84E+05  2.01E+02
 
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
 #CPUT: Total CPU Time in Seconds,       29.861
Stop Time:
Wed Sep 29 14:14:06 CDT 2021
